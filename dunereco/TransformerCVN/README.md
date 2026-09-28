# TransformerCVN

Joint event and prong classifier for DUNE FD: a convolutional encoder turns an
image of the whole event and one image per Pandora track and shower into
embeddings, and a transformer encoder combines them to classify the event and
each prong at once. This directory holds the art integration for the
atmospheric/beam classifier (`art/`: mapper, evaluator and pixel-map producer;
`func/`: data products), introduced in
[PR #128](https://github.com/DUNE/dunereco/pull/128).

`nnbar/` is a separate classification with a similar backbone: the n-nbar
search trains and evaluates its own TransformerCVN network on its own pixel
maps, outside art. Everything for it, from the pixel-map macro to training
and inference, lives in that directory as a standalone tree that the build
does not touch; see `nnbar/README.md`.
