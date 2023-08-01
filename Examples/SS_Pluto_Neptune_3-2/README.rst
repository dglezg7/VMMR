Pluto and Neptune: a 3:2 First Order Resonance
============

--------

Simulation Info:

===================   ============
**Date**              2023/08/01
**Author**            David Graham
**Modules**           SpiNBody
**Approx. runtime**   1 hour (3.8 GHz CPU base clock)
===================   ============

This example produces all potential 1st-3rd order orbital frequency ratios of all planetary pairs in the solar system. I left the directory "Smaller_dEta" for anybody curious to note the convergence of values in relation to the default dEta = 1.0e-2. All bdoies convergence except for the Venus-Mercury, an interesting case for anybody curious to continue investigating.

Expected output
---------------

.. figure:: ResArgPair_Pluto_Neptune.png?raw=True

The Pluto-Neptune 3:2 MMR remains stable up to 1 million years. The Resonant argument on the top panel slowly circulates roughly every 800 kyears while on the bottom panel the angle centers itself at 180 degrees and remains librating.
