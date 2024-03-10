namespace FEM2;

public class Spline
{
    private double h;
    private Point2D[] points;
    private Matrix localMatrix;
    private SparseMatrix globalMatrix;
    private Vector vectorB;
    private Solver slae = default!;
    private Vector w;
    private Vector q;
    private int elementNum;
    private double alpha;
    private double beta;

    public Point2D FirstValue { get => points[0]; } 

    public Spline(int elementNum, string path, (double, double) parametres)
    {

        this.elementNum = elementNum;
        localMatrix = new(4);
        vectorB = new(2 * (elementNum + 1));
        q = new(2 * (elementNum + 1));
        int length;

        (alpha, beta) = parametres;

        using (StreamReader sr = new(path))
        {
            string[] data;
            data = sr.ReadLine().Split(" ");
            length = Convert.ToInt32(data[0]);
            points = new Point2D[length];

            int xIndex = path == "mu" ? 1 : 0, yIndex = (xIndex + 1) % 2;

            for (int i = 0; i < length; i++)
            {
                data = sr.ReadLine().Split(" ").Where(str => str != "").ToArray();
                points[i] = new Point2D(double.Parse(data[xIndex]), double.Parse(data[yIndex]));
            }
        }

        points.OrderBy(point => point.X).ToArray();

        h = (points[^1].X - points[0].X) / elementNum;
        w = new(length);

        w.Fill(1.0);
    }

    public void Compute()
    {
        BuildPortrait();
        AssemblySLAE();
        slae = new(10000, 1e-16);
        slae.SetSLAE(vectorB, globalMatrix);
        slae.Solve();
        Vector.Copy(slae.solution, q);
    }

    private void AssemblySLAE()
    {
        int pointIdx = 0;
        double leftBorder, rightBorder;
        double start = points[0].X;

        for (int i = 0; i < elementNum; i++)
        {
            if (!(pointIdx < points.Length)) break;

            leftBorder = start + i * h;
            rightBorder = start + (i + 1) * h;

            for (;pointIdx < points.Length;)
            {
                if (points[pointIdx].X - rightBorder > 1e-7)
                    break;

                double ksi = getKsi(points[pointIdx].X, i);

                for (int iPsi = 0; iPsi < localMatrix.Size; iPsi++)
                {
                    for (int jPsi = 0; jPsi < localMatrix.Size; jPsi++)
                    {
                        localMatrix[iPsi, jPsi] += w[pointIdx] * getPsi(ksi, iPsi) * getPsi(ksi, jPsi);
                        
                        double FunctionD(double point)
                            => getDPsi(point, iPsi) * getDPsi(point, jPsi);
                        double FunctionDd(double point)
                            => getDdPsi(point, iPsi) * getDdPsi(point, jPsi);

                        localMatrix[iPsi, jPsi] += alpha * Gauss1D(FunctionD, (leftBorder, rightBorder));
                        localMatrix[iPsi, jPsi] += beta * Gauss1D(FunctionDd, (leftBorder, rightBorder));
                    }

                    vectorB[2 * i + iPsi] += w[pointIdx] * getPsi(ksi, iPsi) * points[pointIdx].Y;
                }

                pointIdx++;
            }

            AddLocalMatrixToGlobal(i);
            localMatrix.Clear();
        }
    }

    private void BuildPortrait()
    {
        int size = q.Length;
        globalMatrix = new(size, (size - 2) / 2 * 5 + 1);

        globalMatrix.Ig[0] = 0;
        globalMatrix.Ig[1] = 0;
        globalMatrix.Ig[2] = 1;
        globalMatrix.Ig[3] = 3;
        globalMatrix.Jg[0] = 0;

        for (int i = 4; i < size; i+= 2)
        {
            globalMatrix.Ig[i] = globalMatrix.Ig[i - 1] + 3;
            globalMatrix.Ig[i + 1] = globalMatrix.Ig[i] + 2;
        }
        globalMatrix.Ig[size] = globalMatrix.Ig[size - 1] + 3;

        int col;
        for (int i = 2; i < size; i+= 2)
        {
            col = i - 2;
            for (int j = globalMatrix.Ig[i]; j < globalMatrix.Ig[i + 1]; j++)
            {
                globalMatrix.Jg[j] = col++;
            }

            col = i - 2;
            for (int j = globalMatrix.Ig[i + 1]; j < globalMatrix.Ig[i + 2]; j++)
            {
                globalMatrix.Jg[j] = col++;
            }
        }
    }

    private void AddLocalMatrixToGlobal(int ielem)
    {
        for (int i = 0; i < localMatrix.Size; i++)
            globalMatrix.Di[2 * ielem + i] += localMatrix[i, i];

        for (int i = 1; i < localMatrix.Size; i++)
            for (int j = 0; j < i; j++)
            {
                globalMatrix.Gg[globalMatrix.Ig[2 * ielem + i + 1] - i + j] += localMatrix[i, j];
            }
    }

    public void printSpline(int stepNum)
    {
        using StreamWriter sw = new("..\\..\\..\\resultSpline.txt");
        double point;
        double leftborder;
        double hLoc = h / stepNum;
        double splineValue;
        double start = points[0].X;

        for (int i = 0; i < elementNum; i++)
        {
            leftborder = start + h * i;

            for (int j = 0; j < stepNum; j++)
            {
                splineValue = 0;
                point = leftborder + hLoc * j;

                splineValue = ValueAtPoint(point);

                sw.WriteLine($"{point} {splineValue}");
            }
        }
        splineValue = 0;

        point = points[^1].X;
        splineValue = ValueAtPoint(point);

        sw.WriteLine($"{point} {splineValue}");
    }

    private double getKsi(double point, int elementNum)
        => (point - (points[0].X + h * elementNum)) / h;

    private double getPsi(double ksi, int psiNum)
        => psiNum switch
        {
            0 => 1 - 3 * ksi * ksi + 2 * ksi * ksi * ksi,
            1 => h * (ksi - 2 * ksi * ksi + ksi * ksi * ksi),
            2 => 3 * ksi * ksi - 2 * ksi * ksi * ksi,
            3 => h * (-ksi * ksi + ksi * ksi * ksi),
            _ => throw new NotImplementedException()
        };

    private double getDPsi(double ksi, int psiNum)
        => psiNum switch
        {
            0 => -6 * (ksi - ksi * ksi) / h,
            1 => 1 - 4 * ksi + 3 * ksi * ksi,
            2 => 6 * (ksi - ksi * ksi) / h,
            3 => -2 * ksi + 3 * ksi * ksi,
            _ => throw new NotImplementedException()
        };

    private double getDdPsi(double ksi, int psiNum)
        => psiNum switch
        {
            0 => -6 * (1 - 2 * ksi) / (h * h),
            1 => (-4 + 6 * ksi) / h,
            2 => 6 * (1 - 2 * ksi) / (h * h),
            3 => (-2 + 6 * ksi) / h,
            _ => throw new NotImplementedException()
        };
    
    public static double Gauss1D(Func<double, double> func, (double, double) interval)
    {
        double res = 0;
        double[] p = [-Math.Sqrt(3.0 / 5.0), 0.0, Math.Sqrt(3.0 / 5.0)];
        double[] w = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
        double x0 = interval.Item1;
        double hx = interval.Item2 - interval.Item1;

        for (int i = 0; i < p.Length; i++)
        {
            double point = hx * (p[i] + 1) / 2 + x0;
            res += w[i] * func(point);
        }

        return res * hx / 2.0;
    }

    private int FindElem(double point)
    {
        double start = points[0].X;
        double end = points[^1].X;

        if (start - point > 1e-13)
            return -1;

        if (point - end > 1e-13)
            return -2;

        double rightBorder;

        for (int i = 0; i < elementNum; i++)
        {
            rightBorder = start + h * (i + 1);
            if (rightBorder - point > 1e-14 || rightBorder == point) return i;

        }

        return -1;
    }

    public double ValueAtPoint(double point)
    {
        double res = 0;
        int elem = FindElem(point);
        //if (elem == -1) Console.WriteLine("B < B_first");

        if (elem == -1) return points[0].Y;

        //if (elem == -2) Console.WriteLine("B > B_last");

        if (elem == -2) return points[^1].X / point * (1.0 / points[^1].Y - 1) + 1;

        double ksi = getKsi(point, elem);
        for (int i = 0; i < localMatrix.Size; i++)
        {
            res += q[2 * elem + i] * getPsi(ksi, i);
        }
        return res;
    }

}
