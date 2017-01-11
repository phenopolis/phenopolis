// NUMERIC FORMATTING

function Fixed( s, wid, dec ) {
   // many combinations of possibilities

   // maybe prepare for upcoming truncate
   var z = 1
   if (dec > 0) {
      z /= Math.pow( 10, dec );
      if (s < -z)  s -= 0.5 * z;
      else
         if (s > z)  s += 0.5 * z;
         else
            s = 0;
      }

   // assure a string
   s = "" + s;

   // chop neg, if any
   var neg = 0;
   if (s.charAt(0) == "-") {
      neg = 2;
      s = s.substring( 1, s.length );
      }

   // chop exponent, if any
   var exp = "";
   var e = s.lastIndexOf( "E" );
   if (e < 0)  e = s.lastIndexOf( "e" );
   if (e > -1) {
      exp = s.substring( e, s.length );
      s = s.substring( 0, e );
      }

   // if dec > 0 assure "."; dp == index of "."
   var dp = s.indexOf( ".", 0 );
   if (dp == -1) {
      dp = s.length;
      if (dec > 0) {
         s += ".";
         dp = s.length - 1;
         }
      }

   // assure leading digit
   if (dp == 0) {
      s = '0' + s;
      dp = 1;
      }

   // not enough dec pl?  add 0's
   while ((dec > 0) && ((s.length - dp - 1) < dec))
      s += "0";

   // too many dec pl?  take a substring
   var places = s.length - dp - 1;
   if (places > dec)
      if (dec == 0)
         s = s.substring( 0, dp );
      else
         s = s.substring( 0, dp + dec + 1 );

   // recover exponent, if any
   s += exp;

   // recover neg, if any
   if (neg > 0)
      s = "-" + s;

   // if not enough width, add spaces IN FRONT
   //    too many places?  tough!
   while (s.length < wid)
      s = " " + s;

   return s
   }

function Prb( x ) {
   if (x < 0)  x = 0;
   else
      if (x > 1)  x = 1;
   return x;
   }

function PosV( x ) {
   if (x < 0)  x = -x;
   return x;
   }


// FACTORIALS

function Fact( x ) {
   // x factorial
   var  t=1;
   while (x > 1)
      t *= x--;
   return t;
   }

function LnFact( x ) {
   // ln(x!) by Stirling's formula
   //   see Knuth I: 111
   if (x <= 1)  x = 1;

   if (x < 12)
      return Math.log( Fact(Math.round(x)) );
   else {
      var invx = 1 / x;
      var invx2 = invx * invx;
      var invx3 = invx2 * invx;
      var invx5 = invx3 * invx2;
      var invx7 = invx5 * invx2;

      var sum = ((x + 0.5) * Math.log(x)) - x;
      sum += Math.log(2*Math.PI) / 2;
      sum += (invx / 12) - (invx3 / 360);
      sum += (invx5 / 1260) - (invx7 / 1680);

      return sum;
      }
   }


// COMBINATIONS

function LnComb( n, k ) {
   if ((k == 0) || (k == n))  return 0;
   else
      if ((k > n) || (k < 0))  return -1E38;
      else
         return  (LnFact(n) - LnFact(k) - LnFact(n-k));
   }


// NORMAL

function NormalP( x ) {
   // Abramowitz & Stegun 26.2.19
   var
      d1 = 0.0498673470,
      d2 = 0.0211410061,
      d3 = 0.0032776263,
      d4 = 0.0000380036,
      d5 = 0.0000488906,
      d6 = 0.0000053830,
      a = Math.abs(x),
      t;

   t = 1.0 + a*(d1+a*(d2+a*(d3+a*(d4+a*(d5+a*d6)))));

   // to 16th power
   t *= t;  t *= t;  t *= t;  t *= t;
   t = 1.0 / (t+t);  // the MINUS 16th

   if (x >= 0)  t = 1-t;
   return t;
   }


// BINOMIAL

function g( x ) {
   // Peizer & Pratt 1968, JASA 63: 1416-1456
   var  switchlev = 0.1, z;

   if (x == 0)  z = 1;
   else
      if (x == 1)  z = 0;
      else {

         var d = 1 - x;

         if (Math.abs(d) > switchlev)
            z = (1 - (x * x) + (2 * x * Math.log(x))) / (d * d);
         else {

            z = d / 3;  // first term
            var di = d;  // d**1

            for (var i = 2; i <= 7; i++) {
               di *= d;  // d**i
               z += (2 * di) / ((i+1) * (i+2));
               }
            }
         }
   return z;
   }

function BinomialPF( p, n, k ) {
   // by Normal approximation }
   // Peizer & Pratt 1968, JASA 63: 1416-1456
   var  inv2 = 1/2, inv3 = 1/3, inv6 = 1/6;

   if (k >= n)  z = 1;
   else {
      var q = 1 - p;

      var s = k + inv2;
      var t = n - k - inv2;

      var d1 = s + inv6 - (n + inv3) * p;
      var d2 = q /(s+inv2)  -  p/(t+inv2)  +  (q-inv2)/(n+1);
      d2 = d1 + 0.02 * d2;

      var num = 1  +  q * g(s/(n*p))  +  p * g(t/(n*q));
      var den = (n + inv6) * p * q;
      var z = num / den;
      z = d2 * Math.sqrt(z);
      z = NormalP( z );
      }

   return z;
   }

function BinomTerm( p, n, k ) {
   // for success probability p and n trials
   //     probability of exactly k successes
   return Math.exp( LnComb(n,k)
                    + k * Math.log(p)
                    + (n-k) * Math.log(1-p) );
   }

function BinomialP( p, n, k ) {
   if (n >= 1000)  return BinomialPF( p, n, k );
   else {
      // term-by-term
      if ((k > n) || (p >= 1))  return 1;
      else {
         var  q = 1 - p;
         var  n1p = (n+1) * p;

         var  t = n * Math.log(q);  // k = 0
         var  r = Math.exp(t);
         var  j = 1;
         while (j <= k) {
            t += Math.log( 1 + (n1p - j) / (j * q) );
            r += Math.exp(t);
            j++;
            }

         return  r;
         }
      }
   }

function binomial_cdf(k,n,p) {
   var p = Prb(parseFloat(p));
   var n = PosV(parseInt(n));
   var k = PosV(parseInt(k));
   if (k > n)  k = n;
   var tnk = Fixed(BinomTerm( p, n, k ),8,4);
   var t = BinomialP( p, n, k );
   var pnk = Fixed(t,8,4);
   var qnk = Fixed(1-t,8,4);
   return t;
   }


