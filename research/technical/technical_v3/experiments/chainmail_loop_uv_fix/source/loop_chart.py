"""Closed torus chart. UVs belong to face corners, not welded vertices."""
import numpy as np

def ring_loop_chart(n_major=64, n_minor=16):
    if n_major < 3 or n_minor < 3:
        raise ValueError("A closed tube needs at least three samples per turn")
    i, j = np.meshgrid(np.arange(n_major), np.arange(n_minor), indexing="ij")
    a = np.stack([j/n_minor, i/n_major], axis=-1)
    b = np.stack([j/n_minor, (i+1)/n_major], axis=-1)
    c = np.stack([(j+1)/n_minor, (i+1)/n_major], axis=-1)
    d = np.stack([(j+1)/n_minor, i/n_major], axis=-1)
    return np.stack([a,b,c,d],axis=2).reshape(-1,4,2).astype(np.float32)

