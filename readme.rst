
**************
Raster Retrace
**************

Image tracing utility.


Feature Set
===========

- Curve (re)fitting, using an iterative simplification algorithm: see
  `curve-fit-nd <https://github.com/ideasman42/curve-fit-nd>`__ library.
- Black and white image tracing.
- Corner detection (with angle threshold).
- SVG vector output.

.. note::

   This is an initial release,
   currently this tool works but only loads ``PPM`` images and writes out ``SVG``.

   Support for other image formats is planned.


Examples
========

Examples below use ``TANGENT`` and ``PIXEL`` passes to show the curve fit.

.. figure:: https://cloud.githubusercontent.com/assets/1869379/26520327/6cead016-4313-11e7-9a98-1ec17fdb5a23.png
   :target: https://github.com/ideasman42/raster-retrace-samples/blob/master/output/tauro_2_only_bull.svg

.. figure:: https://cloud.githubusercontent.com/assets/1869379/26520404/42cfb506-4315-11e7-9f76-a83edb73f868.png
   :target: https://github.com/ideasman42/raster-retrace-samples/blob/master/output/tauro_2.svg

.. figure:: https://cloud.githubusercontent.com/assets/1869379/26520321/6049d294-4313-11e7-82a8-9c29e40c3b43.png
   :target: https://github.com/ideasman42/raster-retrace-samples/blob/master/output/jacqueline_face_i.svg

.. figure:: https://cloud.githubusercontent.com/assets/1869379/26520354/1bd0f858-4314-11e7-9f78-604d0fab5f5d.png
   :target: https://github.com/ideasman42/raster-retrace-samples/blob/master/output/blob_simple.svg

.. figure:: https://cloud.githubusercontent.com/assets/1869379/26520322/62e16620-4313-11e7-9a2f-550c015776ee.png
   :target: https://github.com/ideasman42/raster-retrace-samples/blob/master/output/old_guitarist.svg



Usage
=====

.. Output of '--help'

::
   Bitmap image tracing utility

   Options:
       -h, --help   Print help text


   File Options:

       -i, --input FILEPATH   The file path to use for input
       -o, --output FILEPATH  The file path to use for writing


   Tracing Behavior:

       -m, --mode MODE          The method used for tracing the image in [OUTLINE, CENTER], (defaults to OUTLINE).
       -z, --turnpolicy POLICY  Method for extracting outlines [BLACK, WHITE, MAJORITY, MINORITY], (defaults to MAJORITY).


   Curve Evaluation Options:

       Parameters controlling curve evaluation behavior.

       -e, --error PIXELS      The error threshold (defaults to 1.0)
       -t, --simplify PIXELS   Simplify polygon before fitting (defaults to 2.0)
       -c, --corner DEGREES    The corner threshold (`pi` or greater to disable, defaults to 30.0)
       --optimize-exhaustive   When passed, perform exhaustive curve fitting (can be slow!)


   Output Options:

       Generic options for output (format agnostic).

       -s, --scale SCALE    Scale for output, (defaults to 1).
       -p, --passes PASSES  Write extra debug graphics, comma separated list of passes including [PIXEL, PRE_FIT, TANGENT], (defaults to []).
       --pass-scale SCALE   Scale graphic details used in some debug passes, (defaults to 1).


TODO
====

While the basics work, currently there are areas for improvement.

- Support for multiple image formats *(most likely using the piston crate)*.
- Improve bitmap outline extraction method.
- Improve center-line extraction method.
