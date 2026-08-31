"""Optional offline percussion-analysis package.

Import concrete tools from their submodules after installing the ``analysis``
extra. Keeping package import empty prevents SciPy from becoming a dependency
of the native DSP bindings or normal Rack builds.
"""

__all__: list[str] = []
