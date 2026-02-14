int rows = 800; // y
int cols = 1200; // x
float aspectRatio = (float)cols / (float)rows;
int[][] iterMap = new int[rows][cols];

double imagSize = 2;
double realSize = imagSize * aspectRatio;
double imagCenter = 0;
double realCenter = 0;

double realMin = realCenter - 0.5 * realSize;
double realMax = realCenter + 0.5 * realSize;
double imagMin = imagCenter - 0.5 * imagSize;
double imagMax = imagCenter + 0.5 * imagSize;

int maxIter = 0;

double[] realSpace;
double[] imagSpace;

void setup() {
  size(1200, 800);
  colorMode(HSB, 360, 100, 100);
}

void draw() {
  realSpace = linspace(realMin, realMax, cols);
  imagSpace = linspace(imagMax, imagMin, rows);
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      Complex z = fromRect(realSpace[col], imagSpace[row]);
      double prevReal = 0;
      double prevImag = 0;

      for (int k = 0; k < 50; k++) {
        Complex fz = f(z);
        Complex dfz = df(z);
        Complex ddfz = ddf(z);
        Complex denom1 = dfz.power(2).scalarMult(2);
        Complex denom2 = fz.mult(ddfz);
        Complex denominator = denom1.subtract(denom2);
        if (denominator.r < 1e-15) {
          break;
        }
        Complex numerator = fz.mult(dfz).scalarMult(2);
        Complex ratio = numerator.divide(denominator);
        z = z.subtract(ratio);

        double dist = Math.sqrt(Math.pow(z.real() - prevReal, 2) + Math.pow(z.imag() - prevImag, 2));
        if (k > 0 && dist < 1e-06) {
          iterMap[row][col] = k;
          break;
        }
        prevReal = z.real();
        prevImag = z.imag();
      }
    }
  }

  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      if (iterMap[row][col] > maxIter) {
        maxIter = iterMap[row][col];
      }
    }
  }

  int[] histogram = new int[maxIter + 1];
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      histogram[iterMap[row][col]]++;
    }
  }

  int[] cumulative = new int[maxIter + 1];
  cumulative[0] = histogram[0];
  for (int i = 1; i <= maxIter; i++) {
    // cumulative can be normalized by dividing by cumulative[maxIter] (last entry in array)
    cumulative[i] = cumulative[i-1] + histogram[i];
  }
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      int iterations = iterMap[row][col];
      float normalized = (float)cumulative[iterations] / cumulative[maxIter];
      float hue = map(normalized, 0, 1, 62, 310);
      stroke(hue, 90, 76);
      point(col, row); // x, y
    }
  }
  noLoop();
  // println("Framecount: ", frameCount);
  // saveFrame("frame-#####.png");
}

Complex f(Complex z) {
  return z.power(100).subtract(z.power(4)).scalarSum(1);
}

Complex df(Complex z) {
  return z.power(99).scalarMult(100).subtract(z.power(3).scalarMult(4));
}

Complex ddf(Complex z) {
  return z.power(98).scalarMult(9900).subtract(z.power(2).scalarMult(12));
}

// Complex f(Complex z) {
//   return z.power(3).scalarSubtract(1);
// }
//
// Complex df(Complex z) {
//   return z.power(2).scalarMult(3);
// }
//
// Complex ddf(Complex z) {
//   return z.scalarMult(6);
// }

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
    if (other.r < 1e-10) {
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
    if (s < 1e-10) {
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
