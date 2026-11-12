"""Minimal web front end: one form, one report page."""
from flask import Flask, Response, render_template_string, request

from .check import Log, run
from .report import to_html, to_json

PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NMR novelty check</title>
<style>
  body { font: 15px/1.5 system-ui, sans-serif; margin: 2rem auto; max-width: 70rem;
         padding: 0 1rem; color: #1b1f24; }
  h1 { font-size: 1.5rem; } h2 { margin-top: 2rem; font-size: 1.2rem; }
  h3 { font-size: 1.05rem; }
  form { display: grid; grid-template-columns: 11rem 1fr; gap: .5rem .8rem;
         align-items: center; max-width: 48rem; }
  input[type=text] { padding: .4rem .5rem; font: inherit; border: 1px solid #b6bcc4;
                     border-radius: 4px; }
  button { grid-column: 2; justify-self: start; padding: .45rem 1.2rem; font: inherit;
           cursor: pointer; }
  .hint { grid-column: 2; color: #5a636e; font-size: .85rem; margin-top: -.3rem; }
  table { border-collapse: collapse; margin: .5rem 0 1rem; font-size: .88rem;
          display: block; overflow-x: auto; }
  th, td { border: 1px solid #d5d9de; padding: .3rem .55rem; text-align: left;
           vertical-align: top; }
  th { background: #f1f3f5; }
  code { background: #f1f3f5; padding: 0 .25rem; border-radius: 3px;
         word-break: break-all; }
  .error { color: #a4161a; }
  .bar { margin: 1.5rem 0 0; font-size: .9rem; }
</style>
</head>
<body>
<h1>NMR novelty check</h1>
<form method="get" action="/">
  <label for="target">Target or reaction</label>
  <input type="text" id="target" name="target" value="{{ q.target }}"
         placeholder="c1ccccc1" required>
  <div class="hint">SMILES, or reaction SMILES reactants&gt;agents&gt;products</div>
  <label for="reactants">Reactants</label>
  <input type="text" id="reactants" name="reactants" value="{{ q.reactants }}">
  <label for="reagents">Reagents</label>
  <input type="text" id="reagents" name="reagents" value="{{ q.reagents }}">
  <label for="impurities">Possible impurities</label>
  <input type="text" id="impurities" name="impurities" value="{{ q.impurities }}">
  <div class="hint">Several SMILES separated by spaces or commas</div>
  <label for="solvent">NMR solvent</label>
  <input type="text" id="solvent" name="solvent" value="{{ q.solvent }}"
         placeholder="CDCl3">
  <button type="submit">Check</button>
</form>
{% if error %}<p class="error">{{ error }}</p>{% endif %}
{% if body %}
<p class="bar"><a href="{{ json_url }}">Download this report as JSON</a></p>
{{ body|safe }}
{% endif %}
</body>
</html>
"""


def _split(text):
    return [s for s in (text or "").replace(",", " ").split() if s]


def create_app(index, startup_log=None):
    app = Flask(__name__)
    carried = list(startup_log.entries) if startup_log else []

    def do_check(args):
        log = Log()
        log.entries.extend(carried)
        return run(index, args.get("target", "").strip(),
                   reactants=_split(args.get("reactants")),
                   reagents=_split(args.get("reagents")),
                   impurities=_split(args.get("impurities")),
                   solvent=args.get("solvent", "").strip() or None, log=log)

    @app.get("/")
    def home():
        q = {k: request.args.get(k, "") for k in
             ("target", "reactants", "reagents", "impurities", "solvent")}
        body = error = None
        if q["target"].strip():
            try:
                body = to_html(do_check(request.args))
            except ValueError as exc:
                error = str(exc)
        json_url = "/api/check?" + request.query_string.decode()
        return render_template_string(PAGE, q=q, body=body, error=error,
                                      json_url=json_url)

    @app.get("/api/check")
    def api():
        if not request.args.get("target", "").strip():
            return {"error": "target is required"}, 400
        try:
            report = do_check(request.args)
        except ValueError as exc:
            return {"error": str(exc)}, 400
        return Response(to_json(report), mimetype="application/json")

    return app
