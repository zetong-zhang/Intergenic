# Contributing to Intergenic

Contributions are welcome through public GitHub issues and pull requests.

## Before opening an issue

- Reproduce the problem with the latest `main` branch.
- Run `mingw32-make test` or `make test`.
- Include the command line, input format, operating system, compiler version,
  expected result, and observed result.
- Reduce genomic inputs to the smallest reproducible example when possible.

Do not upload confidential, personally identifying, or restricted genomic data.

## Pull requests

1. Open an issue for substantial behavior or file-format changes before starting
   implementation.
2. Keep each pull request focused and preserve unrelated user changes.
3. Add or update automated tests for changed behavior.
4. Run the full test suite on a clean build.
5. Update `README.md`, `CHANGELOG.md`, and user documentation when behavior
   changes.

```powershell
mingw32-make clean all test
```

On Linux or macOS:

```sh
make clean all test
```

Pull requests should explain the user-visible effect, compatibility impact, and
verification steps. New third-party code must have a compatible license and be
documented in `LICENSE`.

## Coding conventions

- C++17, standard library first.
- Keep public behavior in `include/` headers and implementation details in
  `src/`.
- Use the existing coordinate conventions consistently.
- Avoid adding generated binaries, coverage files, or large test outputs to Git.

By contributing, you agree that your contribution is licensed under GPL-3.0-only
and that you will follow [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
