"""Small safeguards shared by crash searches and their independent audits."""


def training_seeds(primary):
    """Always use two different realizations, even for a user-supplied seed."""
    return [primary, 1982 if primary != 1982 else 1944]


def audit_seeds(*primary_seeds):
    """Avoid all supplied fits' training seeds when choosing audit realizations."""
    used = {seed for primary in primary_seeds for seed in training_seeds(primary)}
    choices = (73519, 41273, 69101, 91387, 112337, 130021)
    result = [seed for seed in choices if seed not in used][:2]
    if len(result) != 2:
        raise ValueError("Cannot choose two independent audit seeds")
    return result


def prepare_output(output):
    """Never mix new results with an older trace or overwrite an input snapshot."""
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(
            f"Output is not empty: {output}. Choose a fresh --output directory"
        )
    output.mkdir(parents=True, exist_ok=True)
