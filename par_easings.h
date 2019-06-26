// EASINGS :: https://github.com/prideout/par
// Robert Penner's easing functions.
//
// Distributed under the MIT License, see bottom of file.

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

#ifndef PAR_EASINGS_H
#define PAR_EASINGS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PAR_EASINGS_FLOAT
#define PAR_EASINGS_FLOAT float
#endif

PAR_EASINGS_FLOAT par_easings_linear(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_cubic(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_out_cubic(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_out_cubic(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_quad(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_out_quad(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_out_quad(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_elastic(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_out_elastic(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_out_elastic(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_bounce(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_out_bounce(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_out_bounce(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_back(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_out_back(PAR_EASINGS_FLOAT t);
PAR_EASINGS_FLOAT par_easings_in_out_back(PAR_EASINGS_FLOAT t);

#ifndef PAR_PI
#define PAR_PI (3.14159265359)
#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#define PAR_SQR(a) ((a) * (a))
#endif

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_EASINGS_IMPLEMENTATION

#define PARFLT PAR_EASINGS_FLOAT

PARFLT par_easings__linear(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    return c * t / d + b;
}

PARFLT par_easings__in_cubic(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d;
    return c * t * t * t + b;
}

PARFLT par_easings__out_cubic(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t = t / d - 1;
    return c * (t * t * t + 1) + b;
}

PARFLT par_easings__in_out_cubic(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d / 2;
    if (t < 1) {
        return c / 2 * t * t * t + b;
    }
    t -= 2;
    return c / 2 * (t * t * t + 2) + b;
}

PARFLT par_easings__in_quad(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d;
    return c * t * t + b;
}

PARFLT par_easings__out_quad(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d;
    return -c * t * (t - 2) + b;
}

PARFLT par_easings__in_out_quad(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d / 2;
    if (t < 1) {
        return c / 2 * t * t + b;
    }
    --t;
    return -c / 2 * (t * (t - 2) - 1) + b;
}

PARFLT par_easings__in_quart(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d;
    return c * t * t * t * t + b;
}

PARFLT par_easings__out_quart(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t = t / d - 1;
    return -c * (t * t * t * t - 1) + b;
}

PARFLT par_easings__in_out_quart(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d / 2;
    if (t < 1) {
        return c / 2 * t * t * t * t + b;
    }
    t -= 2;
    return -c / 2 * (t * t * t * t - 2) + b;
}

PARFLT par_easings__in_quint(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d;
    return c * t * t * t * t * t + b;
}

PARFLT par_easings__out_quint(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t = t / d - 1;
    return c * (t * t * t * t * t + 1) + b;
}

PARFLT par_easings__in_out_quint(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d / 2;
    if (t < 1) {
        return c / 2 * t * t * t * t * t + b;
    }
    t -= 2;
    return c / 2 * (t * t * t * t * t + 2) + b;
}

PARFLT par_easings__in_sine(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    return -c * cos(t / d * (PAR_PI / 2)) + c + b;
}

PARFLT par_easings__out_sine(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    return c * sin(t / d * (PAR_PI / 2)) + b;
}

PARFLT par_easings__in_out_sine(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    return -c / 2 * (cos(PAR_PI * t / d) - 1) + b;
}

PARFLT par_easings__in_out_expo(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d / 2;
    if (t < 1) {
        return c / 2 * pow(2, 10 * (t - 1)) + b;
    }
    return c / 2 * (-pow(2, -10 * --t) + 2) + b;
}

PARFLT par_easings__in_circ(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d;
    return -c * (sqrt(1 - t * t) - 1) + b;
}

PARFLT par_easings__out_circ(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t = t / d - 1;
    return c * sqrt(1 - t * t) + b;
}

PARFLT par_easings__in_out_circ(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d / 2;
    if (t < 1) {
        return -c / 2 * (sqrt(1 - t * t) - 1) + b;
    }
    t -= 2;
    return c / 2 * (sqrt(1 - t * t) + 1) + b;
}

PARFLT par_easings__in_elastic(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT a, p, s;
    s = 1.70158;
    p = 0;
    a = c;
    if (!p) {
        p = d * 0.3;
    }
    if (a < fabsf(c)) {
        a = c;
        s = p / 4;
    } else {
        s = p / (2 * PAR_PI) * asin(c / a);
    }
    t -= 1;
    return -(a * pow(2, 10 * t) *
           sin((t * d - s) * (2 * PAR_PI) / p)) + b;
}

PARFLT par_easings__out_elastic(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT a, p, s;
    s = 1.70158;
    p = 0;
    a = c;
    if (!p) {
        p = d * 0.3;
    }
    if (a < fabsf(c)) {
        a = c;
        s = p / 4;
    } else {
        s = p / (2 * PAR_PI) * asin(c / a);
    }
    return a * pow(2, -10 * t) *
           sin((t * d - s) * (2 * PAR_PI) / p) + c + b;
}

PARFLT par_easings__in_out_elastic(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT a, p, s;
    s = 1.70158;
    p = 0;
    a = c;
    if (!p) {
        p = d * (0.3 * 1.5);
    }
    if (a < fabsf(c)) {
        a = c;
        s = p / 4;
    } else {
        s = p / (2 * PAR_PI) * asin(c / a);
    }
    if (t < 1) {
        t -= 1;
        return -0.5 * (a * pow(2, 10 * t) *
               sin((t * d - s) * (2 * PAR_PI) / p)) + b;
    }
    t -= 1;
    return a * pow(2, -10 * t) *
           sin((t * d - s) * (2 * PAR_PI) / p) * 0.5 + c + b;
}

PARFLT par_easings__in_back(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT s = 1.70158;
    t/=d;
    return c*t*t*((s+1)*t - s) + b;
}

PARFLT par_easings__out_back(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT s = 1.70158;
    t=t/d-1;
    return c*(t*t*((s+1)*t + s) + 1) + b;
}

PARFLT par_easings__in_out_back(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT s = 1.70158;
    t/=d/2;
    s*=1.525;
    if (t < 1) { return c/2*(t*t*((s+1)*t - s)) + b; }
    s*=1.525;
    t-=2;
    return (c/2*(t*t*(s+1)*t + s) + 2) + b;
}

PARFLT par_easings__out_bounce(PARFLT t, PARFLT b, PARFLT c, PARFLT d);

PARFLT par_easings__in_bounce(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT v = par_easings__out_bounce(d - t, 0, c, d);
    return c - v + b;
}

PARFLT par_easings__out_bounce(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    t /= d;
    if (t < 1.0 / 2.75) {
        return c * (7.5625 * t * t) + b;
    }
    if (t < 2.0 / 2.75) {
        t -= 1.5 / 2.75;
        return c * (7.5625 * t * t + 0.75) + b;
    }
    if (t < 2.5 / 2.75) {
        t -= 2.25 / 2.75;
        return c * (7.5625 * t * t + 0.9375) + b;
    }
    t -= 2.625 / 2.75;
    return c * (7.5625 * t * t + 0.984375) + b;
}

PARFLT par_easings__in_out_bounce(PARFLT t, PARFLT b, PARFLT c, PARFLT d)
{
    PARFLT v;
    if (t < d / 2) {
        v = par_easings__in_bounce(t * 2, 0, c, d);
        return v * 0.5 + b;
    }
    v = par_easings__out_bounce(t * 2 - d, 0, c, d);
    return v * 0.5 + c * 0.5 + b;
}

PARFLT par_easings_linear(PARFLT t)
{
    return par_easings__linear(t, 0, 1, 1);
}

PARFLT par_easings_in_cubic(PARFLT t)
{
    return par_easings__in_cubic(t, 0, 1, 1);
}

PARFLT par_easings_out_cubic(PARFLT t)
{
    return par_easings__out_cubic(t, 0, 1, 1);
}

PARFLT par_easings_in_out_cubic(PARFLT t)
{
    return par_easings__in_out_cubic(t, 0, 1, 1);
}

PARFLT par_easings_in_quad(PARFLT t)
{
    return par_easings__in_quad(t, 0, 1, 1);
}

PARFLT par_easings_out_quad(PARFLT t)
{
    return par_easings__out_quad(t, 0, 1, 1);
}

PARFLT par_easings_in_out_quad(PARFLT t)
{
    return par_easings__in_out_quad(t, 0, 1, 1);
}

PARFLT par_easings_in_elastic(PARFLT t)
{
    return par_easings__in_elastic(t, 0, 1, 1);
}

PARFLT par_easings_out_elastic(PARFLT t)
{
    return par_easings__out_elastic(t, 0, 1, 1);
}

PARFLT par_easings_in_out_elastic(PARFLT t)
{
    return par_easings__in_out_elastic(t, 0, 1, 1);
}

PARFLT par_easings_in_bounce(PARFLT t)
{
    return par_easings__in_bounce(t, 0, 1, 1);
}

PARFLT par_easings_out_bounce(PARFLT t)
{
    return par_easings__out_bounce(t, 0, 1, 1);
}

PARFLT par_easings_in_out_bounce(PARFLT t)
{
    return par_easings__in_out_bounce(t, 0, 1, 1);
}

PARFLT par_easings_in_back(PARFLT t)
{
    return par_easings__in_back(t, 0, 1, 1);
}

PARFLT par_easings_out_back(PARFLT t)
{
    return par_easings__out_back(t, 0, 1, 1);
}

PARFLT par_easings_in_out_back(PARFLT t)
{
    return par_easings__in_out_back(t, 0, 1, 1);
}

#undef PARFLT

#endif // PAR_EASINGS_IMPLEMENTATION
#endif // PAR_EASINGS_H

// par_easings is distributed under the MIT license:
//
// Copyright (c) 2019 Philip Rideout
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
