CADET
======

.. image:: https://img.shields.io/github/release/modsim/cadet.svg
   :target: https://github.com/modsim/CADET/releases

.. image:: https://github.com/modsim/CADET/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/modsim/CADET/actions/workflows/ci.yml

.. image:: https://anaconda.org/conda-forge/cadet/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/cadet

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8179015.svg
   :target: https://doi.org/10.5281/zenodo.8179015

- **Website (including documentation):** https://cadet.github.io
- **Forum:** https://forum.cadet-web.de
- **Source:** https://github.com/modsim/cadet
- **Bug reports:** https://github.com/modsim/cadet/issues
- **Demo:** https://www.cadet-web.de 
- **Newsletter:** https://cadet-web.de/newsletter/

Installation
------------
CADET can be installed via conda from the ``conda-forge`` channel.

``conda install -c conda-forge cadet``

This requires a working `conda installation <https://docs.anaconda.com/anaconda/install/index.html>`_.

Optionally, use `mamba <https://github.com/mamba-org/mamba>`_ which uses a faster dependency solver than ``conda``.

``mamba install -c conda-forge cadet``

`Additional information <https://cadet.github.io/master/getting_started/installation>`_ and a `tutorial <https://cadet.github.io/master/getting_started/tutorials/breakthrough>`_ are available to guide you through the installation and the first steps of using CADET.

Citing
------------
The development of CADET has been a remarkable collaborative effort, with multiple dedicated individuals contributing their expertise to create a powerful and versatile open-source software tool. Countless hours of hard work and innovation have been invested to provide the scientific community with a valuable resource for chromatography analysis and design. As an open-source project, CADET relies on the support and recognition from users and researchers to thrive. Therefore, we kindly request that any publications, research, or work leveraging the capabilities of CADET acknowledge its creators and their seminal contributions by citing a relevant subset of the following publications.

Cite the software in general:

- **Leweke, S.; von Lieres, E.: Chromatography Analysis and Design Toolkit (CADET)**, Computers and Chemical Engineering 113 (2018), 274–294. `DOI: 10.1016/j.compchemeng.2018.02.025 <https://doi.org/10.1016/j.compchemeng.2018.02.025>`_

Cite numerical algorithms and techniques:

- **von Lieres, E.; Andersson, J.: A fast and accurate solver for the general rate model of column liquid chromatography**, Computers and Chemical Engineering 34,8 (2010), 1180–1191. `DOI: 10.1016/j.compchemeng.2010.03.008 <https://doi.org/10.1016/j.compchemeng.2010.03.008>`_

Cite parameter sensitivities and algorithmic differentiation:

- **Püttmann, A.; Schnittert, S.; Leweke, S.; von Lieres, E.: Utilizing algorithmic differentiation to efficiently compute chromatograms and parameter sensitivities**, Chemical Engineering Science, 139 (2016), 152–162. `DOI: 10.1016/j.ces.2015.08.050 <https://doi.org/10.1016/j.ces.2015.08.050>`_

- **Püttmann, A.; Schnittert, S.; Naumann, U.; von Lieres, E.: Fast and accurate parameter sensitivities for the general rate model of column liquid chromatography**, Computers and Chemical Engineering 56,13 (2013), 46-57. `DOI: 10.1016/j.compchemeng.2013.04.021 <https://doi.org/10.1016/j.compchemeng.2013.04.021>`_

Additionally, we recommend citing the zenodo doi for the specific release that you used in your work to ensure reproducibility.

Ongoing Development
-------------------

We do our best to provide you with a stable API. However, CADET is actively developed and breaking changes can sometimes be unavoidable. For non-developers, it is recommended to upgrade from release to release instead of always working with the most recent commit.

Bugs
----

Please report any bugs that you find `here <https://github.com/modsim/cadet/issues>`_. Or, even better, fork the repository on `GitHub <https://github.com/modsim/cadet>`_ and create a pull request (PR) with the fix. 

Donations
---------

`Donations <https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=FCQ2M89558ZAG>`_ for helping to host, maintain, and further develop the CADET project are highly appreciated.


License
----------

Released under GPL v3. License (see `LICENSE.txt <https://github.com/modsim/CADET/blob/master/LICENSE.txt>`_)::

   Copyright (C) 2004-2023 CADET Authors 
