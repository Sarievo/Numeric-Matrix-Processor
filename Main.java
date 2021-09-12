import java.util.Scanner;

public class Main {
    static Scanner s = new Scanner(System.in);
    public static void main(String[] args) {
        boolean flag = true;
        while (flag) {
            System.out.println("1. Add matrices");
            System.out.println("2. Multiply matrix by a constant");
            System.out.println("3. Multiply matrices");
            System.out.println("4. Transpose matrices");
            System.out.println("5. Calculate a determinant");
            System.out.println("6. Inverse matrix");
            System.out.println("0. Exit");
            System.out.print("Your choice: ");
            int choice = s.nextInt();

            switch (choice) {
                case 1:
                    addMatrices.call(s);
                    break;
                case 2:
                    multiplyMatrix(s);
                    break;
                case 3:
                    multiplyMatrices.call(s);
                    break;
                case 4:
                    transposeMatrices.transpose(s);
                    flag = false;
                    break;
                case 5:
                    calculateADeterminant.call(s);
                    break;
                case 6:
                    inverseMatrix.call(s);
                case 0:
                    flag = false;
                    break;
            }
        }
    }
    private static void printMatrix(double[][] matrix) {
        for (double[] doubles : matrix) {
            for (int i = 0; i < matrix[0].length; i++) {
                System.out.print(doubles[i] + " ");
            }
            System.out.println();
        }
    }
    public static void errorMsg() {
        System.out.println("The operation cannot be performed.");
        System.out.println();
    }
    static class newMatrix {
        private static double[][] single(Scanner s) {
            System.out.print("Enter matrix size: \r");
            int r1 = s.nextInt();
            int c1 = s.nextInt();

            double[][] matrix = new double[r1][c1];
            System.out.println("Enter matrix:");
            for (int i = 0; i < r1; i++) {
                for (int j = 0; j < c1; j++) {
                    matrix[i][j] = Double.parseDouble(s.next());
                }
            }

            return matrix;
        }

        private static double[][] secM1(Scanner s) {
            System.out.print("Enter size of first matrix: \r");
            int r1 = s.nextInt();
            int c1 = s.nextInt();

            double[][] matrix1 = new double[r1][c1];
            System.out.println("Enter first matrix:");
            for (int i = 0; i < r1; i++) {
                for (int j = 0; j < c1; j++) {
                    matrix1[i][j] = Double.parseDouble(s.next());
                }
            }
            return matrix1;
        }

        private static double[][] secM2(Scanner s, int r2, int c2) {
            double[][] matrix2 = new double[r2][c2];
            System.out.println("Enter second matrix:");
            for (int i = 0; i < r2; i++) {
                for (int j = 0; j < c2; j++) {
                    matrix2[i][j] = Double.parseDouble(s.next());
                }
            }
            return matrix2;
        }
    }

    static class addMatrices {
        private static void call(Scanner s) {
            double[][] m1 = newMatrix.secM1(s);

            System.out.print("Enter size of second matrix: \r");
            int r2 = s.nextInt();
            int c2 = s.nextInt();

            if (m1.length == r2 && m1[0].length == c2) {
                double[][] m2 = newMatrix.secM2(s, r2, c2);
                printMatrix(add(m1, m2));
            } else {
                errorMsg();
            }
        }

        private static double[][] add(double[][] m1, double[][] m2) {
            double[][] m = m1.clone();
            for (int i = 0; i < m1.length; i++) {
                for (int j = 0; j < m1[0].length; j++) {
                    m[i][j] += m2[i][j];
                }
            }
            return m;
        }
    }
    private static void multiplyMatrix(Scanner s) {
        double[][] m = newMatrix.single(s);
        System.out.print("Enter constant: \r");
        double c = Double.parseDouble(s.next());
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] *= c;
            }
        }
        printMatrix(m);
    }
    static class multiplyMatrices {
        private static void call(Scanner s) {
            double[][] m1 = newMatrix.secM1(s);

            System.out.print("Enter size of second matrix: \r");
            int r2 = s.nextInt();
            int c2 = s.nextInt();

            if (m1[0].length == r2) {
                double[][] m2 = newMatrix.secM2(s, r2, c2);
                System.out.println("The result is:");
                printMatrix(multiply(m1, m2));
            } else errorMsg();
        }

        private static double[][] multiply(double[][] m1, double[][] m2) {
            double[][] m = new double[m1.length][m2[0].length];
            for (int i = 0; i < m1.length; i++) {
                for (int j = 0; j < m2[0].length; j++) {
                    for (int k = 0; k < m1[0].length; k++) {
                        m[i][j] += m1[i][k] * m2[k][j];
                    }
                }
            }
            return m;
        }
    }
    static class transposeMatrices {
        private static void transpose(Scanner s) {
            boolean flag = true;
            while (flag) {
                System.out.println("1. Main diagonal");
                System.out.println("2. Side diagonal");
                System.out.println("3. Vertical line");
                System.out.println("4. Horizontal line");
                System.out.println("0. Exit");
                System.out.println("Your choice: ");
                int choice = s.nextInt();

                switch (choice) {
                    case 1 -> {
                        double[][] t1 = mainDiagonal(newMatrix.single(s));
                        System.out.println("The result is: ");
                        printMatrix(t1);
                    }
                    case 2 -> {
                        double[][] t2 = sideDiagonal(newMatrix.single(s));
                        System.out.println("The result is: ");
                        printMatrix(t2);
                    }
                    case 3 -> {
                        double[][] t3 = verticalLine(newMatrix.single(s));
                        System.out.println("The result is: ");
                        printMatrix(t3);
                    }
                    case 4 -> {
                        double[][] t4 = horizontalLine(newMatrix.single(s));
                        System.out.println("The result is: ");
                        printMatrix(t4);
                    }
                    case 0 -> flag = false;
                }
            }
        }

        private static double[][] mainDiagonal(double[][] matrix) {
            int row = matrix.length;
            int col = matrix[0].length;

            double[][] newMatrix = new double[col][row];
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    newMatrix[j][i] = matrix[i][j];
                }
            }
            return newMatrix;
        }
        private static double[][] sideDiagonal(double[][] matrix) {
            return horizontalLine(verticalLine(mainDiagonal(matrix)));
        }
        private static double[][] verticalLine(double[][] matrix) {
            int row = matrix.length;
            int col = matrix[0].length;

            double[][] newMatrix = new double[row][col];
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    newMatrix[i][newMatrix[0].length - 1 - j] = matrix[i][j];
                }
            }
            return newMatrix;
        }
        private static double[][] horizontalLine(double[][] matrix) {
            int row = matrix.length;
            int col = matrix[0].length;

            double[][] newMatrix = new double[row][col];
            for (int i = 0; i < row; i++) {
                System.arraycopy(matrix[i], 0, newMatrix[newMatrix.length - 1 - i], 0, col);
            }
            return newMatrix;
        }
    }
    static class calculateADeterminant {
        private static void call(Scanner s) {
            double[][] m = newMatrix.single(s);
            System.out.println("The result is:");
            System.out.println(det(m)); // m in d out
        }
        // return the determinant by the corresponding matrix
        private static double det(double[][] m) {
            if (m.length == 2) { // Trivial
                return m[0][0] * m[1][1] - m[0][1] * m[1][0];
            }

            double sum = 0; // iterate *column times
            for (int i = 0; i < m.length; i++) {
                sum += m[0][i] * cofactor(m, 0, i);
            }
            return sum;
        }

        // return the determinant by the corresponding indices of the corresponding matrix
        private static double cofactor(double[][] m, int i1, int i2) {
            return Math.pow(-1, (i1 + i2) % 2) * det(subMatrix(m, i1, i2));
        }
        // return the subMatrix by the corresponding indices
        private static double[][] subMatrix(double[][] m, int i1, int i2) {
            double[][] sub = new double[m.length - 1][m.length - 1];
            for (int i = 0; i < m.length; i++) {
                for (int j = 0; j < m.length; j++) {
                    if (i < i1) {
                        if (j < i2) {
                            sub[i][j] = m[i][j];
                        } else if (j > i2) {
                            sub[i][j - 1] = m[i][j];
                        }
                    } else if (i > i1) {
                        if (j < i2) {
                            sub[i - 1][j] = m[i][j];
                        } else if (j > i2) {
                            sub[i - 1][j - 1] = m[i][j];
                        }
                    }
                }
            }
            return sub;
        }
    }
    static class inverseMatrix {
        private static void call(Scanner s) {
            double[][] m = newMatrix.single(s);
            double detM = calculateADeterminant.det(m);
            if (detM == 0)  {
                System.out.println("This matrix doesn't have an inverse.");
            } else {
                System.out.println("The result is:");
                printMatrix(inverse(m, detM));
            }
        }
        private static double[][] inverse(double[][] m, double detM) {
            double[][] m1 = new double[m.length][m.length];
            for (int i = 0; i < m.length; i++) { // calculate the adjacent matrix
                for (int j = 0; j < m.length; j++) {
                    m1[i][j] = calculateADeterminant.cofactor(m, i, j);
                }
            }
            m1 = transposeMatrices.mainDiagonal(m1); // main diagonal transpose m1
            for (int i = 0; i < m.length; i++) { // calculate the adjacent matrix
                for (int j = 0; j < m.length; j++) {
                    m1[i][j] = m1[i][j] / detM;
                }
            }
            return m1;
        }
    }
}
