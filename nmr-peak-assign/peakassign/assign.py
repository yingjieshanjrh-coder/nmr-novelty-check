"""Sort measured 1H peaks into product, reactant, solvent and unexplained.

Nothing in here knows about chemistry toolkits: a species is a list of
predicted signals (shift, number of protons, search window) and a measured
peak is a shift with an integral. The steps are

1. fit the product. The scale that turns raw integrals into proton counts is
   whatever makes most product signals come out with the right number of
   protons, so no reference peak has to be picked by hand;
2. take out residual solvent and water;
3. fit each reactant on what is left. A leftover reactant has its own
   concentration, so its signals only have to agree with each other
   (one molar ratio per reactant), not with the product scale;
4. whatever still carries integral is reported, loudly if it is large.

A predicted signal is used at most once, and a peak that is too small for
the signal it would be assigned to is never accepted for it.
"""
from dataclasses import dataclass, field
from functools import lru_cache

MAX_MERGE = 3          # predicted signals that may hide under one multiplet
EXACT_LIMIT = 12       # signals in one crowded region searched exhaustively
EXCESS_WEIGHT = 0.5    # a match that leaves integral unexplained counts half
MIN_SHARE = 0.25       # ...and must explain at least this much of the peak
SINGLET_LIKE = {"", "s", "br", "bs", "brs", "br s", "br. s", "m"}


@dataclass
class Group:
    """One predicted signal."""
    label: str
    shift: float
    protons: int
    window: float
    tight: float = 0.08
    radius: int = 4
    labile: bool = False
    atoms: tuple = ()


@dataclass
class Species:
    key: str
    role: str            # product, reactant or impurity
    name: str
    groups: list
    smiles: str = ""
    atom_map: str = ""


@dataclass
class Settings:
    rel_tol: float = 0.15        # integral tolerance, relative part
    abs_tol: float = 0.10        # ...absolute part, in product protons
    unknown_h: float = 0.5       # unexplained integral that triggers a warning
    residual_alert: float = 0.20  # reactant:product ratio called "a lot"
    min_coverage: float = 0.6    # share of a reactant's signals that must show
    labile_window: float = 1.0   # OH/NH signals wander
    strict_overlap: bool = False


@dataclass
class Part:
    """A share of one peak explained by one species."""
    species: str
    role: str
    name: str
    labels: list
    predicted: list
    protons: float       # nuclei the species puts here; None for plain lines
    h: float             # integral taken, in product protons
    status: str          # ok, excess, split, labile or line
    delta: float
    radius: int = 4
    tight: float = 0.08
    confidence: str = ""


@dataclass
class PeakState:
    index: int
    ppm: float
    integral: float
    multiplicity: str = ""
    h: float = 0.0
    remaining: float = 0.0
    parts: list = field(default_factory=list)
    overlaps: list = field(default_factory=list)   # confirmed species hidden here
    also: list = field(default_factory=list)       # windows that merely cover it
    hints: list = field(default_factory=list)
    conflict: str = ""
    category: str = ""
    left_h: float = 0.0


@dataclass
class SpeciesReport:
    species: Species
    status: str = "absent"       # detected, doubtful, masked or absent
    ratio: float = None          # mol per mol product
    matched: list = field(default_factory=list)
    missing: list = field(default_factory=list)
    masked: list = field(default_factory=list)
    h: float = 0.0
    share: float = 0.0
    alert: bool = False


@dataclass
class Result:
    peaks: list
    unit: float
    scale_source: str
    product: SpeciesReport
    reactants: list
    impurities: list
    hidden: list
    warnings: list
    settings: Settings
    solvent: str = ""

    def category(self, name):
        return [p for p in self.peaks if p.category == name]

    @property
    def unknown_alert(self):
        return bool(self.category("unknown") or self.hidden)


@dataclass
class Fit:
    ratio: float
    matches: list        # (peak index, groups, status, expected raw integral)
    score: float


def _closeness(ppm, group):
    return 1.0 - 0.5 * ((ppm - group.shift) / group.window) ** 2


def _exact(rows):
    """Best choice of one option per peak with no signal used twice."""
    bit = {}
    for row in rows:
        for option in row:
            for g in option[1]:
                bit.setdefault(g, 1 << len(bit))
    masks = [[sum(bit[g] for g in option[1]) for option in row]
             for row in rows]

    @lru_cache(maxsize=None)
    def best(i, used):
        if i == len(rows):
            return 0.0, ()
        score, tail = best(i + 1, used)
        top = (score, (None,) + tail)
        for k, option in enumerate(rows[i]):
            if masks[i][k] & used:
                continue
            score, tail = best(i + 1, used | masks[i][k])
            if option[0] + score > top[0] + 1e-9:
                top = (option[0] + score, (k,) + tail)
        return top

    return best(0, 0)[1]


def _aligned(rows, rank):
    """Same, for crowded regions: peaks and signals are paired in shift
    order only, which keeps the search small."""
    spans = [[(min(rank[g] for g in option[1]), max(rank[g] for g in option[1]))
              for option in row] for row in rows]

    @lru_cache(maxsize=None)
    def best(i, floor):
        if i == len(rows):
            return 0.0, ()
        score, tail = best(i + 1, floor)
        top = (score, (None,) + tail)
        for k, option in enumerate(rows[i]):
            first, last = spans[i][k]
            if first < floor:
                continue
            score, tail = best(i + 1, last + 1)
            if option[0] + score > top[0] + 1e-9:
                top = (option[0] + score, (k,) + tail)
        return top

    return best(0, 0)[1]


def _clusters(rows):
    """Split the peaks into sets that do not compete for any signal."""
    owner, members = {}, {}
    for i, row in enumerate(rows):
        touched = {g for option in row for g in option[1]}
        if not touched:
            continue
        merged = {i}
        for root in {owner[g] for g in touched if g in owner}:
            merged |= members.pop(root)
        members[i] = merged
        for peak in merged:
            for option in rows[peak]:
                for g in option[1]:
                    owner[g] = i
    return [sorted(m) for m in members.values()]


def _solve(avail, groups, near, ratio, noise, cfg):
    rows = []
    for idx, ppm, amount in avail:
        row = []
        for sub in near[idx]:
            expected = ratio * sum(groups[g].protons for g in sub)
            tol = cfg.rel_tol * expected + noise
            if amount < expected - tol:
                continue
            excess = amount > expected + tol
            if excess and expected < MIN_SHARE * amount:
                continue
            gain = sum(groups[g].protons * _closeness(ppm, groups[g])
                       for g in sub)
            if excess:
                gain *= EXCESS_WEIGHT
            row.append((gain, sub, "excess" if excess else "ok", expected))
        rows.append(row)

    by_shift = sorted(range(len(groups)), key=lambda g: -groups[g].shift)
    rank = {g: n for n, g in enumerate(by_shift)}
    score, clean, matches, misfit = 0.0, 0.0, [], 0.0
    for cluster in _clusters(rows):
        cluster.sort(key=lambda i: -avail[i][1])
        local = [rows[i] for i in cluster]
        signals = {g for row in local for option in row for g in option[1]}
        picks = _exact(local) if len(signals) <= EXACT_LIMIT \
            else _aligned(local, rank)
        for i, k in zip(cluster, picks):
            if k is None:
                continue
            gain, sub, status, expected = rows[i][k]
            matches.append((avail[i][0], [groups[g] for g in sub], status,
                            expected))
            score += gain
            if status == "ok":
                clean += gain
                misfit += abs(avail[i][2] - expected) / expected
    return score, clean, matches, misfit


def fit_species(avail, groups, cfg, noise=None, ratio=None):
    """Best one-to-one pairing of a species' signals with the given peaks.

    `avail` holds (peak index, ppm, raw integral still unexplained). With
    `ratio` (raw integral per proton) left open, every ratio suggested by a
    peak/signal pair is tried and the one that explains most protons with
    the right integral wins.
    `noise` is the absolute integral tolerance in raw units; by default it
    scales with the ratio, which is right for the product.
    """
    groups = [g for g in groups if not g.labile]
    near = {}
    for idx, ppm, _ in avail:
        close = sorted((g for g, grp in enumerate(groups)
                        if abs(ppm - grp.shift) <= grp.window),
                       key=lambda g: groups[g].shift)
        # signals merged into one multiplet are neighbours in shift
        near[idx] = [tuple(close[i:i + size])
                     for size in range(1, MAX_MERGE + 1)
                     for i in range(len(close) - size + 1)]

    if ratio is not None:
        candidates = [ratio]
    else:
        raw = sorted(amount / sum(groups[g].protons for g in sub)
                     for idx, _, amount in avail for sub in near[idx])
        candidates = []
        for value in raw:
            if not candidates or value > candidates[-1] * 1.02:
                candidates.append(value)

    def attempt(value):
        tol = cfg.abs_tol * value if noise is None else noise
        score, clean, matches, misfit = _solve(avail, groups, near, value,
                                               tol, cfg)
        # protons explained with the right integral decide; matches that
        # leave area over only break ties
        return (round(clean, 6), round(score, 6), -misfit), \
            Fit(value, matches, score)

    best = None
    for value in candidates:
        trial = attempt(value)
        if trial[1].matches and (best is None or trial[0] > best[0]):
            best = trial
    if best is None:
        return None
    if ratio is None:
        # the winning ratio came from a single peak; average over all clean
        # matches and keep the result if it does not lose anything
        amounts = dict((idx, amount) for idx, _, amount in avail)
        clean = [(amounts[idx], sum(g.protons for g in grps))
                 for idx, grps, status, _ in best[1].matches if status == "ok"]
        if clean:
            refined = attempt(sum(a for a, _ in clean) / sum(n for _, n in clean))
            if refined[0] >= best[0]:
                best = refined
    return best[1]


class _Run:
    def __init__(self, peaks, cfg):
        self.cfg = cfg
        self.states = [PeakState(i, p.ppm, p.integral, p.multiplicity,
                                 remaining=p.integral)
                       for i, p in enumerate(peaks)]
        self.unit = 1.0
        self.used = set()
        self.warnings = []

    @property
    def eps(self):
        return 0.02 * self.unit

    @property
    def noise(self):
        return self.cfg.abs_tol * self.unit

    def avail(self):
        return [(s.index, s.ppm, s.remaining) for s in self.states
                if s.remaining > self.eps]

    def in_window(self, group, window=None):
        window = window or group.window
        return [s for s in self.states if abs(s.ppm - group.shift) <= window]

    def take(self, state, species, groups, status, amount):
        amount = min(amount, state.remaining)
        state.remaining -= amount
        if state.remaining < self.eps:
            amount += state.remaining
            state.remaining = 0.0
        state.parts.append(Part(
            species=species.key, role=species.role, name=species.name,
            labels=[g.label for g in groups],
            predicted=[g.shift for g in groups],
            protons=sum(g.protons for g in groups),
            h=amount / self.unit, status=status,
            delta=max(abs(state.ppm - g.shift) for g in groups),
            radius=min(g.radius for g in groups),
            tight=min(g.tight for g in groups)))
        self.used.update(id(g) for g in groups)

    def apply(self, species, fit):
        for idx, groups, status, expected in fit.matches:
            state = self.states[idx]
            amount = expected if status == "excess" else state.remaining
            self.take(state, species, groups, status, amount)

    # -- product ----------------------------------------------------------

    def product(self, species, reference):
        fixed = None
        if reference:
            near = min(self.states,
                       key=lambda s: abs(s.ppm - reference["ppm"]))
            if abs(near.ppm - reference["ppm"]) > 0.05:
                raise ValueError("no measured peak within 0.05 ppm of the "
                                 f"reference at {reference['ppm']} ppm")
            fixed = near.integral / reference["protons"]
        fit = fit_species(self.avail(), species.groups, self.cfg, ratio=fixed)
        if fixed:
            self.unit, source = fixed, "reference"
        elif fit:
            self.unit, source = fit.ratio, "product"
        else:
            protons = sum(g.protons for g in species.groups) or 1
            self.unit = sum(s.integral for s in self.states) / protons
            source = "total"
            self.warnings.append(
                "没有任何实测峰能与产物的预测信号对应，积分只能按总面积粗略归一化；"
                "请检查产物结构、溶剂，或在 targets.yaml 中用 reference 指定一个参考峰。")
        for state in self.states:
            state.h = state.integral / self.unit
        if fit:
            self.apply(species, fit)
        self.split_signals(species)
        return source

    def split_signals(self, species):
        """A product signal reported as two or three neighbouring peaks."""
        for group in species.groups:
            if group.labile or id(group) in self.used:
                continue
            free = sorted((s for s in self.in_window(group) if not s.parts),
                          key=lambda s: s.ppm)
            expected = self.unit * group.protons
            tol = self.cfg.rel_tol * expected + self.noise
            runs = [free[i:i + size] for size in (2, 3)
                    for i in range(len(free) - size + 1)]
            runs = [r for r in runs
                    if abs(sum(s.remaining for s in r) - expected) <= tol]
            if not runs:
                continue
            run = min(runs, key=lambda r: abs(
                sum(s.ppm * s.remaining for s in r)
                / sum(s.remaining for s in r) - group.shift))
            for state in run:
                self.take(state, species, [group], "split", state.remaining)

    # -- solvent lines ----------------------------------------------------

    def line(self, line):
        cands = [s for s in self.states
                 if abs(s.ppm - line["ppm"]) <= line["window"]]
        if line["kind"] == "water":
            cands = [s for s in cands
                     if s.multiplicity.strip().lower() in SINGLET_LIKE]
        free = [s for s in cands if s.remaining > self.eps]
        if not free:
            for state in cands:
                state.also.append(f"{line['name']}（{line['ppm']:.2f}）")
            return
        state = min(free, key=lambda s: abs(s.ppm - line["ppm"]))
        amount = state.remaining
        state.remaining = 0.0
        state.parts.append(Part(
            species=line["name"], role="solvent", name=line["name"],
            labels=[line["kind"]], predicted=[line["ppm"]], protons=None,
            h=amount / self.unit, status="line",
            delta=abs(state.ppm - line["ppm"])))

    # -- reactants and declared impurities --------------------------------

    def guest(self, species, min_coverage):
        report = SpeciesReport(species)
        fit = fit_species(self.avail(), species.groups, self.cfg,
                          noise=self.noise)
        hit = {id(g) for m in (fit.matches if fit else []) for g in m[1]}
        hosts = {}
        for group in species.groups:
            if group.labile:
                continue
            if id(group) in hit:
                report.matched.append(group)
                continue
            covered = [s for s in self.in_window(group)
                       if s.remaining <= self.eps]
            if covered:
                report.masked.append(group)
                hosts[id(group)] = covered
            else:
                report.missing.append(group)

        seen = len(report.matched) + len(report.missing)
        if not report.matched:
            report.status = "masked" if report.masked and not report.missing \
                else "absent"
        elif len(report.matched) / seen >= min_coverage:
            report.status = "detected"
        else:
            report.status = "doubtful"

        if report.status == "detected":
            report.ratio = fit.ratio / self.unit
            before = sum(s.remaining for s in self.states)
            self.apply(species, fit)
            report.h = (before - sum(s.remaining for s in self.states)) \
                / self.unit
            for group in report.masked:
                host = min(hosts[id(group)],
                           key=lambda s: abs(s.ppm - group.shift))
                host.overlaps.append({
                    "name": species.name, "label": group.label,
                    "shift": group.shift,
                    "h": report.ratio * group.protons})
        elif report.status == "doubtful":
            lost = "、".join(f"{g.label}（{g.shift:.2f}）"
                            for g in report.missing)
            for idx, groups, _, _ in fit.matches:
                self.states[idx].hints.append(
                    f"位移落在 {species.name} 的 {groups[0].label} 预测窗口内，"
                    f"但该物质的其余信号 {lost} 没有出现，未按其归属")
        else:
            for group in report.masked:
                for state in hosts[id(group)]:
                    state.also.append(
                        f"{species.name} {group.label}（{group.shift:.2f}）")
        return report

    # -- exchangeable protons ---------------------------------------------

    def labile(self, species, ratio):
        for group in species.groups:
            if not group.labile:
                continue
            window = max(group.window, self.cfg.labile_window)
            limit = ratio * group.protons * (1 + self.cfg.rel_tol) + self.noise
            cands = [s for s in self.in_window(group, window)
                     if not s.parts and s.remaining <= limit]
            if cands:
                state = min(cands, key=lambda s: abs(s.ppm - group.shift))
                self.take(state, species, [group], "labile", state.remaining)

    # -- verdicts ---------------------------------------------------------

    def classify(self, product):
        open_groups = [g for g in product.groups
                       if not g.labile and id(g) not in self.used]
        hidden = []
        for state in self.states:
            state.left_h = state.remaining / self.unit
            if state.left_h < 0.02:
                state.left_h = 0.0
            species = {p.species for p in state.parts}
            if not state.parts:
                clash = [g for g in open_groups
                         if abs(state.ppm - g.shift) <= g.window]
                if clash:
                    group = min(clash, key=lambda g: abs(state.ppm - g.shift))
                    state.category = "ambiguous"
                    state.conflict = (
                        f"位移符合产物的 {group.label}（预测 {group.shift:.2f} "
                        f"ppm，{group.protons}H），但积分只有 {state.h:.2f}H，"
                        "对不上，未按产物归属")
                elif state.left_h >= self.cfg.unknown_h:
                    state.category = "unknown"
                else:
                    state.category = "minor"
                continue
            if (len(species) > 1 or state.overlaps
                    or (self.cfg.strict_overlap and state.also)):
                state.category = "ambiguous"
            else:
                role = state.parts[0].role
                state.category = "solvent" if role == "impurity" else role
            if state.left_h >= self.cfg.unknown_h:
                hidden.append(state)
        return hidden

    def confidence(self, reports):
        verified = {r.species.key for r in reports
                    if len(r.matched) > 1 and not r.missing}
        for state in self.states:
            for part in state.parts:
                if part.role == "solvent":
                    continue
                sure = (part.status == "ok" and part.radius >= 3
                        and state.category in ("product", "reactant", "solvent")
                        and part.delta <= part.tight)
                if part.role != "product" and part.species not in verified:
                    sure = False
                part.confidence = "高" if sure else "中"


def analyse(peaks, product, reactants, lines=(), impurities=(), cfg=None,
            reference=None, solvent=""):
    """Classify every measured peak. See the module docstring for the steps."""
    cfg = cfg or Settings()
    run = _Run(peaks, cfg)
    source = run.product(product, reference)
    for line in lines:
        run.line(line)
    extra = [run.guest(s, 1.0) for s in impurities]
    guests = [run.guest(s, cfg.min_coverage) for s in reactants]

    run.labile(product, run.unit)
    for report in guests:
        if report.status == "detected":
            run.labile(report.species, report.ratio * run.unit)

    hidden = run.classify(product)
    run.confidence(guests + extra)

    total = sum(s.h for s in run.states)
    main = SpeciesReport(product, status="detected", ratio=1.0)
    for group in product.groups:
        if id(group) in run.used:
            main.matched.append(group)
        elif not group.labile:
            main.missing.append(group)
    reports = [main] + guests + extra
    for report in reports:
        report.h = sum(p.h for s in run.states for p in s.parts
                       if p.species == report.species.key)
        report.share = report.h / total if total else 0.0
    if not main.matched:
        main.status = "absent"
    for report in guests:
        report.alert = (report.status == "detected"
                        and report.ratio >= cfg.residual_alert)

    return Result(peaks=run.states, unit=run.unit, scale_source=source,
                  product=main, reactants=guests, impurities=extra,
                  hidden=hidden, warnings=run.warnings, settings=cfg,
                  solvent=solvent)
