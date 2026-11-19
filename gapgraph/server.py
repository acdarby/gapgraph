"""Tiny JSON + static server for the editor UI (stdlib only)."""
import json
import os
import traceback
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
from threading import Lock

STATIC = os.path.join(os.path.dirname(__file__), "static")


def make_handler(session, export_prefix):
    lock = Lock()

    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, body, ctype="application/json"):
            data = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            u = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            try:
                if u.path in ("/", "/index.html"):
                    with open(os.path.join(STATIC, "index.html"), "rb") as fh:
                        return self._send(200, fh.read(), "text/html; charset=utf-8")
                if u.path == "/api/overview":
                    with lock:
                        return self._send(200, session.overview())
                if u.path == "/api/region":
                    return self._send(200, session.region(q["contig"], int(q["start"]), int(q["end"])))
                if u.path == "/api/end":
                    return self._send(200, session.end_detail(q["contig"], q["side"]))
                if u.path == "/api/hic":
                    with lock:
                        return self._send(200, session.hic_matrix() or {})
                if u.path == "/api/problems":
                    if "contig" in q:
                        return self._send(200, session.problems(q["contig"], q["type"], q.get("thr")))
                    return self._send(200, session.problems_summary(q["type"], q.get("thr")))
                if u.path == "/api/next":
                    return self._send(200, session.next_problem(q["contig"], int(q["pos"]), q["type"],
                                                                q.get("thr"), int(q.get("dir", 1))) or {})
                if u.path == "/api/templates":
                    return self._send(200, session.templates(q["contig"], int(q["start"]), int(q["end"])))
                if u.path == "/api/joinview":
                    req = json.loads(q["req"])
                    return self._send(200, session.joinview(req, int(q["start"]), int(q["end"])))
                if u.path == "/api/joinalign":
                    return self._send(200, session.joinalign(json.loads(q["req"])))
                if u.path == "/api/coverage":
                    out = {}
                    for plat, st in session.stores.items():
                        out[plat] = st.coverage(q["contig"], session.lengths[q["contig"]], 300)
                    return self._send(200, out)
                self._send(404, {"error": "not found"})
            except Exception as e:
                traceback.print_exc()
                self._send(400, {"error": str(e)})

        def do_POST(self):
            u = urlparse(self.path)
            n = int(self.headers.get("Content-Length", 0))
            req = json.loads(self.rfile.read(n) or b"{}")
            try:
                with lock:
                    if u.path == "/api/op":
                        res = session.op(req)
                        return self._send(200, {"ok": True, "result": res})
                    if u.path == "/api/export":
                        return self._send(200, {"ok": True, "files": session.project.export(export_prefix)})
                self._send(404, {"error": "not found"})
            except Exception as e:
                self._send(400, {"ok": False, "error": str(e)})

    return H


def serve(session, host="127.0.0.1", port=8765, export_prefix="gapgraph_out"):
    httpd = ThreadingHTTPServer((host, port), make_handler(session, export_prefix))
    print(f"gapgraph editor on http://{host}:{port}  (Ctrl-C to stop)")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
