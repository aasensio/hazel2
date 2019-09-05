.. include:: ../hazel_name
.. _regularization:


Regularization
==============

Sometimes it is useful to introduce some regularization during the inversion process. This is
added as additional terms to the :math:`\chi^2` merit function. For the moment, only the following
regularizations are possible.

l2-value
--------

This type of regularization penalizes values of the parameters that get far from a certain user-defined
value. It is defined by adding the following in the ``[[[Regularization]]]`` section of any atmosphere and
any parameter, passing a list with three elements:

::

    [[[Regularization]]]
    v      = l2-value, 1.0, 5.0

The first element of the list is the one defining the regularization type. The second defines
the value of the regularization hyperparameter, while the last element defines the value that we
want to use in the regularization. A term like the following is added to the merit function:

.. math::

    \chi^2_\mathrm{modified} = \chi^2 + \lambda |v-x|^2

where in the previous case :math:`\lambda=1` and :math:`x=5`. You can see the effect of the regularization
by adding a very large value of :math:`\lambda` and checking that the code shifts the parameter to this
specific value.

l2-gradient
-----------

Not yet implemented