
**************
Raster Retrace
**************

Image tracing utility.


.. note::

   This is an initial release,
   currently this tool works but only loads ``PPM`` images and writes out ``SVG``.

   Support for other image formats is planned.


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

       -m, --mode MODE          The method used for tracing the image in [OUTLINE, CENTER], (defaults to CENTER).
       -z, --turnpolicy POLICY  Method for extracting outlines
                                [BLACK, WHITE, MAJORITY, MINORITY], (defaults to MAJORITY).


   Curve Evaluation Options:

       Parameters controlling curve evaluation behavior.

       -e, --error PIXELS      The error threshold (defaults to 1.0)
       -t, --simplify PIXELS   Simplify polygon before fitting (defaults to 2.0)
       -c, --corner DEGREES    The corner threshold (`pi` or greater to disable, defaults to 30.0)
       --optimize-exhaustive   When passed, perform exhaustive curve fitting (can be slow!)


   Output Options:

       Generic options for output (format agnostic).

       -s, --scale SCALE    Scale for output, (defaults to 1).
       -p, --passes PASSES  Write extra debug passes,
                            comma separated list of passes including [PIXEL, PRE_FIT], (defaults to []).

TODO
====

While the basics work, currently there are areas for improvement.

- Support for multiple image formats *(most likely using the piston crate)*.
- Improve bitmap outline extraction method.
- Improve center-line extraction method.
