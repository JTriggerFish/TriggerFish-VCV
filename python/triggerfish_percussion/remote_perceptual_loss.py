"""Optional SSH scoring transport; the actual synthesizer stays local C++/Wasm."""

import base64
import json
import subprocess

import numpy as np
from .perceptual_fit_losses import ScalarAudioLoss


class RemoteJtfsLoss(ScalarAudioLoss):
    def __init__(self, reference, rate, host, directory, python):
        # Paths come from our explicit development command, never sample tags.
        if any(c in directory + python for c in "\n\r'\""):
            raise ValueError("Remote paths must not contain shell quoting")
        command = (
            f"PYTHONPATH='{directory}/deps:{directory}/python' "
            f"'{python}' '{directory}/perceptual_loss_worker.py'"
        )
        self.process = subprocess.Popen(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host, command],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf8",
        )
        try:
            self.specification = self.request("initialize", reference, rate=rate)[
                "specification"
            ]
        except Exception:
            self.close()
            raise

    def request(self, command, samples, **values):
        pcm = base64.b64encode(np.asarray(samples, dtype="<f4").tobytes()).decode(
            "ascii"
        )
        self.process.stdin.write(
            json.dumps(dict(command=command, pcm=pcm, **values)) + "\n"
        )
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("JTFS worker exited without a response")
        result = json.loads(line)
        if "error" in result:
            raise RuntimeError(result["error"])
        return result

    def score(self, samples):
        return self.request("score", samples)["score"]

    def close(self):
        self.process.stdin.close()
        try:
            self.process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            self.process.terminate()
            self.process.wait(timeout=5)
        self.process.stdout.close()
