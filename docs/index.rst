.. quentem_mechenex documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: qm5.png


Welcome to qm5's documentation!
=========================================================
qm5 is a package that will solve approximatelly for the following Hamiltonian:

.. math:: \hat{H} = E_{ion} + \sum_{p,q} \sum_{\sigma}  h_{p,q} \hat{a}^{\dagger} \hat{a}_{q,\sigma} + \frac{1}{2} \sum_{p,q,r,s} \sum_{\sigma} V_{p,q,r,s} \hat{a}_{p,\sigma}^{\dagger} \hat{a}_{r,\sigma} \hat{a}_{s,\sigma} \hat{a}_{q,\sigma} 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api
   Noble_Gas_Model
   HartreeFock
   HartreeFockFast


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
