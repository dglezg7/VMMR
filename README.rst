VMMR Calculator: Virtual Mean Motion Resonance Calculator
============

Plots Resonant Arguments of Simulations currently exclusive to VPLanet: https://github.com/VirtualPlanetaryLaboratory/vplanet

Run the following terminal commands to simulate any of the following posted examples or your own work.
First, run the simulator (vplanet) to produce results.
second, run vmmr.py. It takes the results to detect Resonant Arguments currently up to 3rd order.

.. code-block:: bash

    # 1. Run the VPLanet simulation
    vplanet vpl.in


.. code-block:: bash

    # 2. Run the calculator to produce the png files.
    python vmmr.py png
    # If you want to observe the next image, press "x" on the currently observed image.
    # To produce a pdf
    python vmmr.py pdf