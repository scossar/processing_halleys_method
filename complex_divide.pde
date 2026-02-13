int rows = 1000;
int cols = rows;
// float scale = 20;
// float scaleAdjust;
int[][] iterMap = new int[rows][cols];
// Complex[][] Z = new Complex[rows][cols];

double realMin = -1;
double realMax = 1;
double imagMin = -1;
double imagMax = 1;

int maxIter = 0;
int[] cumulative;

void setup() {
  size(1000, 1000);
  colorMode(HSB, 360, 100, 100);
  // scaleAdjust = width/scale;
  double[] x = linspace(realMin, realMax, width);
  double[] y = linspace(imagMax, imagMin, height);

  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      Complex z = fromRect(x[i], y[j]);
      double gPrev = 0; // hack

      for (int k = 0; k < 100; k++) {
        Complex fz = f(z);
        Complex dfz = df(z);
        Complex ddfz = ddf(z);
        Complex denom1 = dfz.power(2).scalarMult(2);
        Complex denom2 = fz.mult(ddfz);
        Complex denominator = denom1.subtract(denom2);
        if (denominator.r < 1e-15) {
          break;
        }
        Complex numerator2 = fz.mult(dfz).scalarMult(2);
        Complex numerator = z.subtract(numerator2);
        z = numerator.divide(denominator);

        double g = Math.pow(z.r, 2);
        if (k > 0 && Math.abs(g - gPrev) < 1e-06) {
          iterMap[i][j] = k;
          break;
        }
        gPrev = g;
      }
    }
  }

  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      if (iterMap[i][j] > maxIter) {
        maxIter = iterMap[i][j];
      }
    }
  }

  int[] histogram = new int[maxIter + 1];
  for (int i = 0; i < iterMap.length; i++) {
    for (int j = 0; j < iterMap[i].length; j++) {
      histogram[iterMap[i][j]]++;
    }
  }
  cumulative = new int[maxIter + 1];
  cumulative[0] = histogram[0];
  for (int i = 1; i <= maxIter; i++) {
    // cumulative can be normalized by dividing by cumulative[maxIter] (last entry in array)
    cumulative[i] = cumulative[i-1] + histogram[i];
  }
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      int iter = iterMap[i][j];
      float normalized = (float)cumulative[iter] / cumulative[maxIter];
      float hue = map(normalized, 0, 1, 40, 360);
      stroke(hue, 100, 100);
      point(i, j);
    }
  }
}

void draw() {
  noLoop();
  // translate(width/2, height/2);
  // scale(1, -1);
  // background(217, 2, 94);
  // stroke(217, 90, 19);
  // strokeWeight(4);
  // point(0, 0);
}

Complex f(Complex z) {
  return z.power(3).scalarSubtract(1);
}

Complex df(Complex z) {
  return z.power(2).scalarMult(3);
}

Complex ddf(Complex z) {
  return z.scalarMult(6);
}

double[] linspace(double start, double end, int num) {
  double[] result = new double[num];
  if (num == 1) {
    result[0] = start;
    return result;
  }

  double step = (end - start) / (num - 1);
  for (int i = 0; i < num; i++) {
    result[i] = start + i * step;
  }
  return result;
}

Complex fromRect(double real, double imag) {
  double r = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
  double theta = Math.atan2(imag, real);
  return new Complex(r, theta);
}

class Complex {
  double r, theta;

  Complex(double r, double theta) {
    this.r = r;
    this.theta = theta;
  }


  Complex mult(Complex other) {
    return new Complex(r * other.r, theta + other.theta);
  }

  Complex divide(Complex other) {
    if (other.r < 1e-10) {  // guessing here, not sure about this.
      return new Complex(Double.POSITIVE_INFINITY, Double.NaN);
    }
    return new Complex (r / other.r, theta - other.theta);
  }

  Complex sum(Complex other) {
    double real = this.real() + other.real();
    double imag = this.imag() + other.imag();
    double r = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
    double theta = Math.atan2(imag, real);
    return new Complex(r, theta);
  }

  Complex subtract(Complex other) {
    double real = this.real() - other.real();
    double imag = this.imag() - other.imag();
    double r = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
    double theta = Math.atan2(imag, real);
    return new Complex(r, theta);
  }

  Complex power(double n) {
    return new Complex(Math.pow(r, n), theta * n);
  }

  Complex scalarSum(double s) {
    double real = this.real() + s;
    double imag = this.imag();
    double r = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
    double theta = Math.atan2(imag, real);
    return new Complex(r, theta);
  }

  Complex scalarSubtract(double s) {
    double real = this.real() - s;
    double imag = this.imag();
    double r = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
    double theta = Math.atan2(imag, real);
    return new Complex(r, theta);
  }

  Complex scalarMult(double s) {
    return new Complex(r * s, theta);
  }

  Complex scalarDivide(double s) {
    if (s < 1e-10) {  // guessing here, not sure about this.
      return new Complex(Double.POSITIVE_INFINITY, theta);
    }
    return new Complex(r / s, theta);
  }

  double real() {
    return r * Math.cos(theta);
  }
  double imag() {
    return r * Math.sin(theta);
  }
}
