"""Exact UI renderer transport; no Python reimplementation of recipe macros."""

import base64
import json
import subprocess
from pathlib import Path

import numpy as np


class WorkbenchRenderer:
    def __init__(self, node: str, target: str, root: Path):
        self.process = subprocess.Popen(
            [node, str(root / "tools/workbench_fit_bridge.mjs")],
            cwd=root,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )
        try:
            self.metadata = self.request(command="initialize", id=target)
        except Exception:
            self.close()
            raise
        self.reference = self.decode(self.metadata.pop("pcm"))
        self.sample_rate = self.metadata["reference"]["sampleRate"]
        self.initial = dict(
            zip(
                [item["key"] for item in self.metadata["descriptors"]],
                self.metadata["values"],
            )
        )

    @staticmethod
    def decode(pcm):
        return np.frombuffer(base64.b64decode(pcm), dtype="<f4").astype(np.float64)

    def request(self, **request):
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("Workbench renderer exited without a response")
        response = json.loads(line)
        if "error" in response:
            raise RuntimeError(response["error"])
        return response

    def render(self, parameters, seconds, seed=None):
        request = dict(parameters=parameters, seconds=seconds)
        if seed is not None:
            request["seed"] = seed
        return self.decode(self.request(**request)["pcm"])

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=10)
        self.process.stdout.close()
