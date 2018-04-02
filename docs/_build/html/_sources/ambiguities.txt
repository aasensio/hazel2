Ambiguities in the Hanle effect in the saturation regime
========================================================

In the saturation regime of the Hanle effect, Stokes :math:`Q` and
:math:`U` are insensitive to the field strength, but are sensitive to
the geometry of the field. For a :math:`J=0 \to J=1` transition, the
linear polarization can be written as:

.. math::

   \begin{aligned}
   Q &=& \frac{q}{2} \left( 3 \cos^2 \theta_B-1 \right) \sin^2\Theta_B \cos 2\Phi_B \nonumber \\
   U &=& \frac{q}{2} \left( 3 \cos^2 \theta_B-1 \right) \sin^2\Theta_B \sin 2\Phi_B.\end{aligned}

These expressions contain a mixture of angles to make it clear that the
polarization amplitude depends on both the angle between the vertical
and the magnetic field and between the magnetic field and the
line-of-sight (LOS).

The coordinates of the magnetic field vector :math:`\mathbf{B}` in the
reference system of the vertical and the reference system of the LOS
are:

.. math::

   \begin{aligned}
   \mathbf{B} &=& B \left(\sin \theta_B \cos \phi_B \mathbf{i}+\sin \theta_B \sin \phi_B \mathbf{j}+\cos \theta_B \mathbf{k} \right) \nonumber \\
   \mathbf{B} &=& B \left(\sin \Theta_B \cos \Phi_B \mathbf{i}'+\sin \Theta_B \sin \Phi_B \mathbf{j}'+\cos \Theta_B \mathbf{k}' \right),\end{aligned}

where the unit vectors are related by a simple rotation:

.. math::

   \begin{aligned}
   \mathbf{i}' &=& \cos \theta \mathbf{i} - \sin \theta \mathbf{k} \nonumber \\
   \mathbf{k}' &=& \sin \theta \mathbf{i} + \cos \theta \mathbf{k}.\end{aligned}

Introducing these relations on the expression for the magnetic field,
we find that the following has to be fulfilled, given that the magnetic
field vector is the same in both reference systems:

.. math::

   \begin{aligned}
   \sin \theta_B \cos \phi_B &=& \sin \Theta_B \cos \Phi_B \cos \theta + \cos \Theta_B + \sin \theta \nonumber \\
   \sin \theta_B \sin \phi_B &=& \sin \Theta_B \sin \Phi_B \nonumber \\
   \cos \theta_B &=& \cos \Theta_B \cos \theta - \sin \Theta_B \cos \Phi_B \sin \theta.\end{aligned}

Solving the previous three equations in the two directions, we find the
following transformations between the angles in the vertical reference
system and the LOS reference system:

.. math::

   \begin{aligned}
   \cos \Theta_B &=& \cos\theta \cos\theta_B + \sin\theta \sin\theta_B \cos\phi_B \nonumber \\
   \sin \Theta_B &=& +\sqrt{1-\cos^2\Theta_B} \nonumber \\
   \cos \Phi_B &=& \frac{\cos\theta \sin\theta_B \cos\phi_B - \cos\theta_B \sin\theta}{\sin \Theta_B} \nonumber \\
   \sin \Phi_B &=& \frac{\sin\theta_B \sin\phi_B}{\sin\Theta_B}\end{aligned}

 and

.. math::

   \begin{aligned}
   \cos \theta_B &=& \cos\theta \cos\Theta_B - \sin\theta \sin\Theta_B \cos\Phi_B \nonumber \\
   \sin \theta_B &=& +\sqrt{1-\cos^2\theta_B} \nonumber \\
   \cos \phi_B &=& \frac{\cos\theta \sin\Theta_B \cos\Phi_B + \cos\Theta_B \sin\theta}{\sin \theta_B} \nonumber \\
   \sin \phi_B &=& \frac{\sin\Theta_B \sin\Phi_B}{\sin\theta_B}.\end{aligned}

Note that, since :math:`\Theta_B \in [0,\pi]`, we can safely use the
square root and take the positive value. In order to transform from one
reference system to the other, we can compute the inclination easily by
inverting the sinus or the cosinus. However, the situation is different
for the azimuth, because the range of variation is :math:`[-\pi,\pi]`.
Therefore, one has to compute the cosinus and the sinus separately and
the decide which is the correct quadrant fo the angle in terms of the
signs of both quantities.

Four possible kinds of ambiguities can exist for the Stokes :math:`Q`
and :math:`U` parameters. The idea is that :math:`\Phi_B` can be
modified and still obtain the same :math:`Q` and :math:`U` by properly
adjusting the value of :math:`\Theta_B`. It is clear that, given that
the term that can be used to compensate for the change in the azimuth on
the LOS reference system is the same for Stokes :math:`Q` and :math:`U`,
we can only compensate for changes in the sign. Therefore, we have the
following potential ambiguities:

.. math::

   \begin{aligned}
   \Phi_B' &=& \Phi_B \nonumber \\
   \Phi_B' &=& \Phi_B -\pi/2 \nonumber \\
   \Phi_B' &=& \Phi_B + \pi/2 \nonumber \\
   \Phi_B' &=& \Phi_B + \pi.\end{aligned}

For each case, we have to compute the value of :math:`\Theta_B'` that
keeps the value of :math:`Q` and :math:`U` unchanged. Therefore, once we
find a solution to the inversion problem in the form of the pair
:math:`(\theta_B,\phi_B)`, we can find the remaining solutions in the
saturation regime following the recipes that we present now. Remember
that, unless one knows the polarity of the field, or in other words, the
sign :math:`\cos\Theta_B`, the number of potential ambiguous solutions
is 8. If the polarity of the field is known, the number is typically
reduced to 4 (or 2 if no 90\ :math:`^\circ` ambiguity is present).

:math:`\Phi_B’=\Phi_B`
----------------------

Under this change, we have that

.. math:: \cos 2\Phi_B' = \cos 2\Phi_B, \quad \sin 2\Phi_B' = \sin 2\Phi_B, \quad \cos \Phi_B' = \cos \Phi_B, \quad \sin \Phi_B' = \sin \Phi_B.

 Making use of the previous relations between the angles wrt to the
vertical and the LOS, we have to solve the following equation:

.. math:: \left( 3 \cos^2\theta_B'-1 \right) \sin^2 \Theta_B' = \left( 3 \cos^2\theta_B-1 \right) \sin^2 \Theta_B,

 which can be written as:

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta - \sin\theta \sin\Theta_B' \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

After some algebra and doing the substitution :math:`t=\sin\Theta_B'`,
we end up with the following equation to be solved:

.. math:: A t^4 + Bt^2 + C t^3 \sqrt{1-t^2} = K,

where

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \cos^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& -6 \cos\theta \sin\theta \cos \Phi_B \nonumber \\
   K &=& \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

The previous equation can be solved if we make the change of variables
:math:`t=\pm \sqrt{Z}`, resulting in:

.. math:: (C^2+A^2) Z^4 + (-C^2+2AB) Z^3 + (-2AK+B^2) Z^2 - 2BKZ + K^2 = 0.

This polynomial of 4-th order can have four different solutions. From
these solutions, we have to take only the real solutions which are
larger than 0, given the range of variation of :math:`\Theta_B`:

.. math:: t \in \mathbb{R}, \qquad 0 \leq t \leq 1.

Once the solutions for :math:`t` are found, we make
:math:`\Theta_B' = \arcsin t`. Note that, for a fixed value of
:math:`t`, two values of :math:`\Theta_B'` are possible. We choose the
correct one by evaluating the expressions for :math:`Q` and :math:`U`
and testing which of the two possible choices give the values equal (or
very similar) to the original ones.

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B)` to the vertical reference
system.

:math:`\Phi_B’=\Phi_B+\pi`
--------------------------


Under this change, we have:

.. math:: \cos 2\Phi_B' = \cos 2\Phi_B, \quad \sin 2\Phi_B' = \sin 2\Phi_B, \quad \cos \Phi_B' = -\cos \Phi_B, \quad \sin \Phi_B' = -\sin \Phi_B.

Following the same approach, we have to solve for :math:`\Theta_B'` in

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta + \sin\theta \sin\Theta_B' \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

The solution are obtained as the roots of the same equations as before
but now

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \cos^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& 6 \cos\theta \sin\theta \cos \Phi_B \nonumber \\
   K &=& \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B+\pi)` to the vertical
reference system.

:math:`\Phi_B’=\Phi_B+\pi/2`
----------------------------

Under this change, we have:

.. math:: \cos 2\Phi_B' = -\cos 2\Phi_B, \quad \sin 2\Phi_B' = -\sin 2\Phi_B, \quad \cos \Phi_B' = -\sin \Phi_B, \quad \sin \Phi_B' = \cos \Phi_B.

Following the same approach, we have to solve for :math:`\Theta_B'` in

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta + \sin\theta \sin\Theta_B' \sin\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

The solution are obtained as the roots of the same equations as before
but now

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \sin^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& 6 \cos\theta \sin\theta \sin \Phi_B \nonumber \\
   K &=& -\left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B+\pi/2)` to the vertical
reference system.

:math:`\Phi_B’=\Phi_B-\pi/2`
----------------------------

Under this change, we have:

.. math:: \cos 2\Phi_B' = -\cos 2\Phi_B, \quad \sin 2\Phi_B' = -\sin 2\Phi_B, \quad \cos \Phi_B' = \sin \Phi_B, \quad \sin \Phi_B' = -\cos \Phi_B.

Following the same approach, we have to solve for :math:`\Theta_B'` in

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta + \sin\theta \sin\Theta_B' \sin\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

The solution are obtained as the roots of the same equations as before
but now

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \sin^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& -6 \cos\theta \sin\theta \sin \Phi_B \nonumber \\
   K &=& -\left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B-\pi/2)` to the vertical
reference system.
