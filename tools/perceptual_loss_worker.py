"""Persistent JTFS scoring worker. Audio is transient input, never written here."""

import base64
import json
import sys

import numpy as np
import torch
from triggerfish_percussion.perceptual_fit_losses import JtfsLoss


def decode(message):
    return np.frombuffer(base64.b64decode(message["pcm"]), dtype="<f4")


def main():
    torch.set_num_threads(2)
    loss = None
    for line in sys.stdin:
        try:
            message = json.loads(line)
            if message["command"] == "initialize":
                loss = JtfsLoss(decode(message), message["rate"], device="cuda")
                result = dict(specification=loss.specification)
            elif message["command"] == "score" and loss is not None:
                result = dict(score=loss.score(decode(message)))
            else:
                raise ValueError("Initialize before scoring")
        except Exception as error:
            result = dict(error=f"{type(error).__name__}: {error}")
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
