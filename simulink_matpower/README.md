SimulinkMATPOWER - A Simulink to [MATPOWER][1] interface
===================================================

This toolbox called [SimulinkMATPOWER][2] enables you to use
[MATPOWER][1] in [Simulink][3]. Furthermore, it adds a control layer for
tap changers if this feature is activated.

The idea is the following: MATPOWER provides us with a power flow
solver. Simulink provides a graphic interface. Let's bring these two
together and enable people to use MATPOWER in Simulink.

>   **Note:** SimulinkMATPOWER is part of [MATPOWER Extras][8],
    which is included in MATPOWER zip distributions in the
    `extras/simulink_matpower` directory.

How it works
------------

What this toolbox does is that it enables you to manipulate a MATPOWER
casefile (`mpc`) from Simulink, run a power flow solver on the `mpc`
from Simulink and then extract the power flow solution in Simulink.
Given that the toolbox purely interconnects MATPOWER with Simulink
everything is based on the standard MATPOWER casefile. We write values
into the `mpc` from Simulink, we run the power flow on the `mpc` from
Simulink and read values from the `mpc` to Simulink.

A Simulink Library was created for this toolbox. The library is called
`SimulinkMATPOWERbase.slx` and can be found in the root folder of the
library. In there you can find the block that enables using MATPOWER's
power flow solver in Simulink (the block "ACPF"). There are also other
blocks that help you interconnect your own blocks with MATPOWER (e.g.
blocks that get voltage magnitudes etc).


Setup
-----

1. [Install MATPOWER][4].
2. Add the toolbox directory **with sub-directories** to the Matlab path.

   > _**Note:** Step 2 is necessary **only** if SimulinkMATPOWER
      was not already included with your MATPOWER distribution (in
      `extras/simulink_matpower`)._

3. Go into `ac_testbed_setup.m` and change the code in line 16 to load
   the MATPOWER case file you want to work with.
4. Run `ac_testbed_setup` (in `example.slx` it is automatically run by
   Simulink as an install_function).
5. Run the example `example.slx` to see whether it works and to get a
   first idea of how this toolbox works.


How to Use the Toolbox
----------------------

1. Drag and drop or copy the ACPF block from the library to your Simulink
   model. This block calls the MATPOWER power flow solver from Simulink.

2. Writing data into the MATPOWER casefile from Simulink. The ACPF block
   has inputs. These inputs correspond to the structure of
   a `mpc` file. For example,

   - `active_power_generation_in` is written into `mpc.gen(:,PG)`
   - `voltage_magnitude_in` is written into `mpc.bus(:,VM)`

   You can use these inputs to parametrise the `mpc` as you want.
   Everything that you do not want to actively change will still need to
   be given to the ACPF block every time the power flow is solved.
   
   Check the example how you can write these values into the ACPF block
   easily (take the output of the ACPF block, store it by one time step
   and then send that to the inputs).

   For example, this enbales you to simulate a closed-loop control
   system in which a controller changes reactive power setpoints of
   generators. Then the power flow is solved to model the power grid and
   the solution of the power flow is given to the controller which then
   again updates the setpoints and so on.

3. Running the MATPOWER power flow solver on the `mpc`.

   This is done automatically at every timestep when you run the Simulink
   model.

   _(optional)_ In line 30 of `ac_testbed_setup.m` you can enable tap
   changer control by setting `mpm.oltc.enable=1`. The default is 0. When
   enabled the following happens: The power flow is solved using
   MATPOWER, then for all generators with tap changers it is checked
   whether or not the secondary voltage is outside the deadband. When
   the voltage is outside the deadband the tap ratio is changed (which
   is exactly what a tap changer would do) and the power flow is solved
   once more. Again all secondary voltages are checked and if needed the
   tap changer position are changed again. This repeats until not more
   tap changes occur. Overall, this calculates the steady-state power
   flow solution including tap changer behaviour.

4. Getting data from the solved power flow into Simulink (aka reading
   values from the `mpc` file after the power flow was solved).

   For this there are several different blocks in the library. For
   example, blocks that extract the voltage magnitudes, the power flows,
   the losses, etc. Drag and drop or copy them into your Simulink file.
   Check the example to see where they should be located (in the
   triggered subsystem "Extracting Power Flow Results").

Checkout the `example.slx` to get your project going.


Important Information
---------------------

1. DEFINE TAP CHANGERS. The file `find_OLTC_indexes` creates a vector
   `idx_branches_in_oltc` with the indexes of the branches that are
   transformers with on-load-tap-changers (OLTCs). Adjust and
   parametrize this file to make it find the indexes of the OLTCs in
   your grid model.

2. If `idx_branches_in_oltc` is empty, then the Simulink subsystem
   extracting the tap ratios must be commented out.

3. The MATPOWER case file `mpc` is inside the `mpm` structure (`mpm.mpc`).


[Citing SimulinkMATPOWER][5]
----------------------------

We request that publications derived from the use of the SimulinkMATPOWER
toolbox explicitly acknowledge that fact by citing the following paper
which initiated the development of this toolbox.

>   Lukas Ortmann, Jean Maeght, Patrick Panciatici, Florian Döfler,
    Saverio Bolognani, "Online Feedback Optimization for Transmission
    Grid Operation," 2022.
    doi: [10.48550/arXiv.2212.07795][6]


Questions
---------

If you have any questions or are facing difficulties, then please reach
out to Lukas Ortmann at lukas.ortmann@web.de.


License
-------

SimulinkMATPOWER is distributed under the [3-clause BSD license][7].

----
[1]: https://matpower.org
[2]: https://github.com/Lukas738/SimulinkMATPOWER
[3]: https://www.mathworks.com/products/simulink.html
[4]: https://matpower.org/about/get-started/
[5]: CITATION
[6]: https://doi.org/10.48550/arXiv.2212.07795
[7]: LICENSE
[8]: https://github.com/MATPOWER/matpower-extras
