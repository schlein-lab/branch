# LLmap dependency pin

BRANCH consumes the LLmap *classical* alignment engine (gap-affine WFA2/Gotoh
pairwise alignment + seed-chain-extend mapping) in place of the former
minimap2/ksw2 dependency. See `src/align/llmap_align_backend.{hpp,cpp}` for the
single seam.

## Pinned commit

- Repository: https://github.com/schlein-lab/LLmap
- Branch: `feat/optional-llm-checkpoint`
- Commit: `0a84749a5923b06df9e111883f0ac58b3b8080f2`

This commit makes LLmap's LLM-consult (claude_agent) path optional via
`LLMAP_ENABLE_CLAUDE`, which BRANCH requires so it can link `llmap_classical`
without the agent/HTTP/CUDA-sandbox stack.

## Building the LLmap libraries BRANCH links against

Configure LLmap with the AI/agent paths off (BRANCH needs only the classical
alignment libraries — `llmap_classical`, `llmap_annot`, `llmap_checkpoint`,
`llmap_core`):

```sh
cmake -S <llmap> -B <llmap>/build-min \
  -DCMAKE_BUILD_TYPE=Release \
  -DLLMAP_ENABLE_CUDA=OFF -DLLMAP_ENABLE_FOUNDATION=OFF \
  -DLLMAP_ENABLE_FAISS=OFF -DLLMAP_ENABLE_CLAUDE=OFF \
  -DLLMAP_ENABLE_TESTS=OFF -DLLMAP_ENABLE_BENCH=OFF \
  -DLLMAP_USE_NATIVE_ARCH=OFF
cmake --build <llmap>/build-min --target llmap_classical -j8
```

## Pointing BRANCH at LLmap

```sh
cmake -S . -B build \
  -DBRANCH_LLMAP_DIR=<llmap> \
  -DBRANCH_LLMAP_BUILD_DIR=<llmap>/build-min
```

`BRANCH_LLMAP_DIR` supplies headers (expects `<dir>/src`); `BRANCH_LLMAP_BUILD_DIR`
supplies the static libraries (expects `<dir>/src` or `<dir>/lib`). The build is
fully offline — there is no network FetchContent.
