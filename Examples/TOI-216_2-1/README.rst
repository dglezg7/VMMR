TOI-216 (2 Confirmed Planets): a 2:1 First Order Resonance
============

--------

Simulation Info:

===================   ============
**Date**              2023/08/01
**Author**            David Graham
**Modules**           SpiNBody
**Approx. runtime**   13 seconds (3.8 GHz CPU base clock)
===================   ============

This example produces the 2:1 MMR of planets b and c in the system. The initial values of TOI-216 were obtained from this publication (Kipping 2019): https://arxiv.org/pdf/1902.03900.pdf 
The paper produces detailed information about the orbital elements. I calculated the Mean Anomaly knowing the planets' transit times and other orbital elements.

To run this example
-------------------

.. code-block:: bash
    # Run vplanet to produce the results
    vplanet vpl.in
    
    # Run the calculator to produce the png files.
    python vmmr.py png
    
    # If you want to observe the next image, press "x" on the currently observed image.

    # To produce a pdf
    python vmmr.py pdf

Expected output
---------------

.. figure:: ResArgPair_c_b.png?raw=True

The two planets show a case of high-amplitude libration on the top panel for up to 1 million years.
