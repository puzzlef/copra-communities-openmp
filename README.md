Multi-threaded OpenMP-based Community OVerlap PRopagation Algorithm ([COPRA]) for
[community detection].

This is an implementation of a label-propagation based community detection
algorithm called **Community OVerlap PRopagation Algorithm (COPRA)**. Unlike
[RAK], this algorithm uses *multiple labels per vertex*, with each label having an
associated *belonging coefficient* (which sums to `1`). The algorithm is as follows.

1. Each vertex initializes as its own community (belonging=1).
2. Each iteration, a vertex collects labels from its neighborhood.
3. I excludes a vertex's own labels, although not explicitly mentioned in paper.
4. The collected labels are scaled by edge weights for weighted graph.
5. Each vertex picks labels above a certain threshold.
6. This threshold is inversely proportional to the max. number of labels.
7. If all labels are below threshold, pick a random best label.
8. I make a vertex join its own community if it has no labels (not mentioned).
9. Selected labels are normalized such that belonging coefficient sums to 1.
10. Repeat from 2 until convergence.

The authors of this paper mention to use a *mimimum vertices per community count*
to detect convergence. However i do not find it to be helpful. I instead use a
similar convergence condition as *RAK*, that is **count the number of vertices that**
**change thier best label**. Once this count falls below a certain fraction
(tolerance), i consider the algorithm to have converged. The authors also
mention using a *sychrounous* version of the algorithm, where labels of each
vertex is dependent only upon labels in previous iteration. However, i find
**asynchronous** approach to be converge faster (labels of each vertex can be
dependent upon labels in current iteration). As we focus of finding disjoint
communities, i consider the **best label of each vertex** as the final result.

I continue with *OpenMP implementation* of **COPRA** algorithm for community
detection. Each thread is given a *separate hashtable*, which it can use for
choosing a set of labels from its neighbors above a certain threshold (by
weight). If no label is above threshold, we pick only the best weighted label,
and if not we simply join our own community. Hashtables are *allocated separately*
for better performance (instead of storing them contiguously on a vector).
OpenMP schedule is `auto` now, we can later choose the best if we need.

[![](https://i.imgur.com/sSbAM0p.png)][sheetp]

<br>
<br>

Similar to [previous experiment], i vary the **tolerance** from `0.1` to `0.0001`,
and adjust the **max.** **number of labels** from `1` to `32`.

[![](https://i.imgur.com/Kd2bsQS.png)][sheetp]

[![](https://i.imgur.com/yWsaYUv.png)][sheetp]

[![](https://i.imgur.com/ywYY5tR.png)][sheetp]

<br>
<br>

On average, *OpenMP approaches are faster* than sequential, and also seem to
*achieve better modularity* than sequential. We also observe that *using a single*
*label per vertex seems to yield both best performance and best modularity* on
average. This is particulary true in case of road networks, but *not so in case*
*of other classes of graphs* (web graphs, social networks, collaboration
networks). For example, *web graphs* such as `web-Stanford` and `web-BerkStan` achieve
best modularity with **max. labels** of `8`, `web-Google` and `web-NotreDame` does best
with **max. labels** of `16`/`32`. **Max. labels** of `4`-`16` would be a good choice
for such graphs.

[![](https://i.imgur.com/NZ2JJIW.png)][sheetp]

<br>
<br>

In addition it seems that on average, *making the tolerance tighter than 0.01*
*has no beneficial effect on modularity*. However, *tighter tolerance does not*
*help with* graphs such as `coAuthorsDBLP` and *social networks*. It seems a
**tolerance** of `0.01` would be a good choice on average.

[![](https://i.imgur.com/9CraoiB.png)][sheetp]

<br>
<br>

Both **RAK** and **COPRA** approaches can obtain *disconnected communities*. This issue
can be resolved by splitting such communities into separate communities in a
*post-processing step*. I do **not** do that here. We may do it when comparing these
approaches.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].


[COPRA]: https://arxiv.org/abs/0910.5516
[RAK]: https://arxiv.org/abs/0709.2938
[community detection]: https://en.wikipedia.org/wiki/Community_search
[previous experiment]: https://github.com/puzzlef/copra-communities-seq
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu

<br>

```bash
$ g++ -std=c++17 -O3 main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 3985272 [directed] {} (symmetricize)
# OMP_NUM_THREADS=12
# [-0.000497 modularity] noop
# [00204.560 ms; 0003 iters.; 0.836576819 modularity] copraSeqStatic {labels=01, tolerance=1e-01}
# [00047.392 ms; 0003 iters.; 0.840981662 modularity] copraOmpStatic {labels=01, tolerance=1e-01}
# [00244.446 ms; 0003 iters.; 0.838712633 modularity] copraSeqStatic {labels=02, tolerance=1e-01}
# [00048.273 ms; 0003 iters.; 0.846463919 modularity] copraOmpStatic {labels=02, tolerance=1e-01}
# [00378.092 ms; 0003 iters.; 0.875163972 modularity] copraSeqStatic {labels=04, tolerance=1e-01}
# [00065.594 ms; 0003 iters.; 0.877589226 modularity] copraOmpStatic {labels=04, tolerance=1e-01}
# [00618.397 ms; 0003 iters.; 0.878541648 modularity] copraSeqStatic {labels=08, tolerance=1e-01}
# [00098.211 ms; 0003 iters.; 0.881777525 modularity] copraOmpStatic {labels=08, tolerance=1e-01}
# [00641.359 ms; 0003 iters.; 0.881162941 modularity] copraSeqStatic {labels=16, tolerance=1e-01}
# [00109.541 ms; 0003 iters.; 0.876401544 modularity] copraOmpStatic {labels=16, tolerance=1e-01}
# [01145.591 ms; 0003 iters.; 0.866800904 modularity] copraSeqStatic {labels=32, tolerance=1e-01}
# [00188.139 ms; 0003 iters.; 0.860131264 modularity] copraOmpStatic {labels=32, tolerance=1e-01}
# [00204.643 ms; 0003 iters.; 0.836576819 modularity] copraSeqStatic {labels=01, tolerance=5e-02}
# [00046.144 ms; 0003 iters.; 0.845682204 modularity] copraOmpStatic {labels=01, tolerance=5e-02}
# [00245.197 ms; 0003 iters.; 0.838712633 modularity] copraSeqStatic {labels=02, tolerance=5e-02}
# [00047.560 ms; 0003 iters.; 0.846338987 modularity] copraOmpStatic {labels=02, tolerance=5e-02}
# [00472.478 ms; 0004 iters.; 0.887068212 modularity] copraSeqStatic {labels=04, tolerance=5e-02}
# [00078.316 ms; 0004 iters.; 0.885929525 modularity] copraOmpStatic {labels=04, tolerance=5e-02}
# [00946.169 ms; 0005 iters.; 0.892077446 modularity] copraSeqStatic {labels=08, tolerance=5e-02}
# [00120.187 ms; 0004 iters.; 0.885639489 modularity] copraOmpStatic {labels=08, tolerance=5e-02}
# [00951.614 ms; 0005 iters.; 0.897701085 modularity] copraSeqStatic {labels=16, tolerance=5e-02}
# [00134.461 ms; 0004 iters.; 0.882408142 modularity] copraOmpStatic {labels=16, tolerance=5e-02}
# [01723.266 ms; 0005 iters.; 0.885300398 modularity] copraSeqStatic {labels=32, tolerance=5e-02}
# ...
```

<br>
<br>


## References

- [Finding overlapping communities in networks by label propagation; Steve Gregory (2010)](https://iopscience.iop.org/article/10.1088/1367-2630/12/10/103018)
- [Near linear time algorithm to detect community structures in large-scale networks; Usha Nandini Raghavan et al. (2007)](https://arxiv.org/abs/0709.2938)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [How to import VSCode keybindings into Visual Studio?](https://stackoverflow.com/a/62417446/1413259)
- [Configure X11 Forwarding with PuTTY and Xming](https://www.centlinux.com/2019/01/configure-x11-forwarding-putty-xming-windows.html)
- [Installing snap on CentOS](https://snapcraft.io/docs/installing-snap-on-centos)

<br>
<br>


[![](https://i.imgur.com/7GLy9tb.jpg)](https://www.youtube.com/watch?v=L-ZBWLYGSuY)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/562704096.svg)](https://zenodo.org/badge/latestdoi/562704096)


[gist]: https://gist.github.com/wolfram77/fcc480d4d549c05cb5482f3cf838efdc
[charts]: https://imgur.com/a/7FQbvW9
[sheets]: https://docs.google.com/spreadsheets/d/1LHb5bFGDATB9NrY1QEt1rHLBu9BGf6b0N3Z44f_s6xQ/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vRmWbXRz-A5z4oSPY_vSy6sZUze-ZT0z79IflsvSNqe7CPbW_EaLHScbawY6DcLE-_fdtqOTmaZE1KI/pubhtml
