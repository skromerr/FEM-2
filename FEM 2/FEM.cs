using System.Drawing;

namespace UMFCourseProject;

public class FEM
{
    public delegate double Basis(Point2D point);

    private const double mu0 = 4 * Math.PI * 1e-7;
    private Grid grid;
    private SparseMatrix globalMatrix;
    private Vector globalVector;
    private Solver slae;
    private Matrix alphas;
    private Vector localVector;
    private Matrix stiffnessMatrix;
    private Matrix massMatrix;
    private Point2D[] vertices;
    private Vector solution;
    private FirstCondition[] firstConditions;
    private SecondCondition[] secondConditions;
    private Basis[] basis;
    private Basis[] dbasisdx;
    private Basis[] dbasisdy;
    private string condPath => "conditions.txt";

    public FEM(Grid grid)
    {
        this.grid = grid;
        alphas = new(3);

        // квадратичный
        //basis = [Psi1, Psi2, Psi3, Psi4, Psi5, Psi6];
        //dbasisdx = [dPsi1dx, dPsi2dx, dPsi3dx, dPsi4dx, dPsi5dx, dPsi6dx];
        //dbasisdy = [dPsi1dy, dPsi2dy, dPsi3dy, dPsi4dy, dPsi5dy, dPsi6dy];

        // линейный
        basis = [_l1, _l2, _l3];
        dbasisdx = [_dl1dx, _dl2dx, _dl3dx];
        dbasisdy = [_dl1dy, _dl2dy, _dl3dy];


        stiffnessMatrix = new(basis.Length);
        localVector = new(basis.Length);
        slae = new Solver(1000, 1e-16);



        vertices = new Point2D[3];
        solution = new Vector(grid.Nodes.Count);

        SetBoundaryConditions();
    }

    private void SetBoundaryConditions()
    {
        firstConditions = new FirstCondition[grid.FirstConditionNodes.Count];
        int i = 0;
        foreach (var node in grid.FirstConditionNodes)
            firstConditions[i++] = new(grid.Nodes[node], node);
    }

    public void Compute()
    {

        BuildPortrait();
        AssemblyGlobalMatrix();
        AccountFirstConditionsWithBigNumber();

        slae.SetSLAE(globalVector, globalMatrix);
        slae.CGM();
        Vector.Copy(slae.solution, solution);

        PrintSolution();

    }


    public void AccountFirstConditions()
    {
        foreach (var fc in firstConditions)
        {
            globalMatrix.Di[fc.NodeNumber] = 1;
            globalVector[fc.NodeNumber] = 0;
            for (int i = globalMatrix.Ig[fc.NodeNumber]; i < globalMatrix.Ig[fc.NodeNumber + 1]; i++)
            {
                globalMatrix.Gg[i] = 0;
            }
            for (int i = fc.NodeNumber + 1; i < globalMatrix.Size; i++)
            {
                for (int j = globalMatrix.Ig[i]; j < globalMatrix.Ig[i + 1]; j++)
                {
                    if (globalMatrix.Jg[j] == fc.NodeNumber)
                    {
                        globalMatrix.Gg[j] = 0;
                    }
                }
            }
        }
    }

    public void AccountFirstConditionsWithBigNumber()
    {
        foreach (var fc in firstConditions)
        {
            globalMatrix.Di[fc.NodeNumber] = 1e30;
            globalVector[fc.NodeNumber] = 0;
        }
    }

    public void BuildPortrait()
    {
        var list = new HashSet<int>[grid.Nodes.Count].Select(_ => new HashSet<int>()).ToList();
        foreach (var element in grid.Elements)
        {
            foreach (var pos in element)
            {
                foreach (var node in element)
                {
                    if (pos > node)
                    {
                        list[pos].Add(node);
                    }
                }
            }
        }

        int count = list.Sum(childList => childList.Count);

        globalMatrix = new(grid.Nodes.Count, count);
        globalVector = new(grid.Nodes.Count);

        globalMatrix.Ig[0] = 0;

        for (int i = 0; i < list.Count; i++)
            globalMatrix.Ig[i + 1] = globalMatrix.Ig[i] + list[i].Count;

        int k = 0;

        foreach (var childList in list)
        {
            foreach (var value in childList.Order())
            {
                globalMatrix.Jg[k++] = value;
            }
        }
    }

    private void AddElement(int i, int j, double value)
    {
        if (i == j)
        {
            globalMatrix.Di[i] += value;
            return;
        }

        for (int icol = globalMatrix.Ig[i]; icol < globalMatrix.Ig[i + 1]; icol++)
        {
            if (globalMatrix.Jg[icol] == j)
            {
                globalMatrix.Gg[icol] += value;
                return;
            }
        }
    }

    private void AssemblyGlobalMatrix()
    {
        globalMatrix.Clear();
        globalVector.Fill(0);

        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
            vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
            vertices[2] = grid.Nodes[grid.Elements[ielem][2]];

            AssemblyLocalMatrixes();

            stiffnessMatrix = 1 / (grid.Materials[grid.Elements[ielem][^1]].Mu) * stiffnessMatrix;

            for (int i = 0; i < basis.Length; i++)
                for (int j = 0; j < basis.Length; j++)
                    AddElement(grid.Elements[ielem][i], grid.Elements[ielem][j], stiffnessMatrix[i, j]);

            AssemblyLocallVector(ielem);

            AddElementToVector(ielem);

            stiffnessMatrix.Clear();
            localVector.Fill(0);
        }
    }

    void AssemblyLocallVector(int ielem)
    {
        double jValue = grid.Materials[grid.Elements[ielem][^1]].J;

        for (int i = 0; i < basis.Length; i++)
        {
            localVector[i] = GaussTriangle(basis[i].Invoke);
        }

        localVector = jValue / 10.0 * Math.Abs(DeterminantD()) * localVector;
    }

    private void AddElementToVector(int ielem)
    {
        for (int i = 0; i < basis.Length; i++)
        {
            globalVector[grid.Elements[ielem][i]] += localVector[i];
        }
    }

    public double GaussTriangle(Func<Point2D, double> func)
    {
        const double x1a = 0.873821971016996;
        const double x1b = 0.063089014491502;
        const double x2a = 0.501426509658179;
        const double x2b = 0.249286745170910;
        const double x3a = 0.636502499121399;
        const double x3b = 0.310352451033785;
        const double x3c = 0.053145049844816;
        const double w1 = 0.050844906370207;
        const double w2 = 0.116786275726379;
        const double w3 = 0.082851075618374;
        double[] p1 = { x1a, x1b, x1b, x2a, x2b, x2b, x3a, x3b, x3a, x3c, x3b, x3c };
        double[] p2 = { x1b, x1a, x1b, x2b, x2a, x2b, x3b, x3a, x3c, x3a, x3c, x3b };
        double[] w = { w1, w1, w1, w2, w2, w2, w3, w3, w3, w3, w3, w3 };
        double res = 0;

        for (int i = 0; i < w.Length; i++)
        {
            Point2D point = new();
            point.X = (1 - p1[i] - p2[i]) * vertices[0].X + p1[i] * vertices[1].X + p2[i] * vertices[2].X;
            point.Y = (1 - p1[i] - p2[i]) * vertices[0].Y + p1[i] * vertices[1].Y + p2[i] * vertices[2].Y;

            res += func(point) * w[i] * 0.5;
        }

        return res;
    }

    private double DeterminantD()
         => (vertices[1].X - vertices[0].X) * (vertices[2].Y - vertices[0].Y) -
            (vertices[2].X - vertices[0].X) * (vertices[1].Y - vertices[0].Y);

    private void CalcuclateAlphas()
    {
        double dD = DeterminantD();

        alphas[0, 0] = (vertices[1].X * vertices[2].Y - vertices[2].X * vertices[1].Y) / dD;
        alphas[0, 1] = (vertices[1].Y - vertices[2].Y) / dD;
        alphas[0, 2] = (vertices[2].X - vertices[1].X) / dD;

        alphas[1, 0] = (vertices[2].X * vertices[0].Y - vertices[0].X * vertices[2].Y) / dD;
        alphas[1, 1] = (vertices[2].Y - vertices[0].Y) / dD;
        alphas[1, 2] = (vertices[0].X - vertices[2].X) / dD;

        alphas[2, 0] = (vertices[0].X * vertices[1].Y - vertices[1].X * vertices[0].Y) / dD;
        alphas[2, 1] = (vertices[0].Y - vertices[1].Y) / dD;
        alphas[2, 2] = (vertices[1].X - vertices[0].X) / dD;
    }

    private void AssemblyLocalMatrixes()
    {
        double dD = Math.Abs(DeterminantD());
        CalcuclateAlphas();

        for (int i = 0; i < stiffnessMatrix.Size; i++)
            for (int j = 0; j < stiffnessMatrix.Size; j++)
            {
                Func<Point2D, double> func = Sum(Mult(dbasisdx[i], dbasisdx[j]), Mult(dbasisdy[i], dbasisdy[j])).Invoke;
                stiffnessMatrix[i, j] = GaussTriangle(func);
            }
        stiffnessMatrix = dD * stiffnessMatrix;
    }

    private (double, double, double) getL(Point2D point)
    {
        double l1 = alphas[0, 0] + alphas[0, 1] * point.X + alphas[0, 2] * point.Y;
        double l2 = alphas[1, 0] + alphas[1, 1] * point.X + alphas[1, 2] * point.Y;
        double l3 = alphas[2, 0] + alphas[2, 1] * point.X + alphas[2, 2] * point.Y;

        return (l1, l2, l3);
    }

    private double basisss(Point2D point, int numPsi)
    {
        (var l1, var l2, var l3) = getL(point);

        switch (numPsi)
        {
            case 0:
                return l1 * (2 * l1 - 1);
            case 1:
                return l2 * (2 * l2 - 1);
            case 2:
                return l3 * (2 * l3 - 1);
            case 3:
                return 4 * l1 * l2;
            case 4:
                return 4 * l2 * l3;
            case 5:
                return 4 * l3 * l1;
            default:
                return 0;
        }
    }

    private double Psi1(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);

        return l1 * (2 * l1 - 1);
    }

    private double Psi2(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);

        return l2 * (2 * l2 - 1);
    }

    private double Psi3(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return l3 * (2 * l3 - 1);
    }

    private double Psi4(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * l1 * l2;
    }
    private double Psi5(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * l2 * l3;
    }

    private double Psi6(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * l3 * l1;
    }

    private double dPsi1dx(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);

        return 4 * alphas[0, 1] * l1 - alphas[0, 1];
    }

    private double dPsi2dx(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);

        return 4 * alphas[1, 1] * l2 - alphas[1, 1];
    }

    private double dPsi3dx(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * alphas[2, 1] * l3 - alphas[2, 1];
    }

    private double dPsi4dx(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * (alphas[0, 1] * l2 + alphas[1, 1] * l1);
    }
    private double dPsi5dx(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * (alphas[1, 1] * l3 + alphas[2, 1] * l2);
    }

    private double dPsi6dx(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * (alphas[2, 1] * l1 + alphas[0, 1] * l3);
    }

    private double dPsi1dy(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);

        return 4 * alphas[0, 2] * l1 - alphas[0, 2];
    }

    private double dPsi2dy(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);

        return 4 * alphas[1, 2] * l2 - alphas[1, 2];
    }

    private double dPsi3dy(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * alphas[2, 2] * l3 - alphas[2, 2];
    }

    private double dPsi4dy(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * (alphas[0, 2] * l2 + alphas[1, 2] * l1);
    }
    private double dPsi5dy(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * (alphas[1, 2] * l3 + alphas[2, 2] * l2);
    }

    private double dPsi6dy(Point2D point)
    {
        (var l1, var l2, var l3) = getL(point);
        return 4 * (alphas[2, 2] * l1 + alphas[0, 2] * l3);
    }

    private double _l1(Point2D point)
        => getL(point).Item1;

    private double _l2(Point2D point)
    => getL(point).Item2;

    private double _l3(Point2D point)
    => getL(point).Item3;

    private double _dl1dx(Point2D point)
    => alphas[0, 1];

    private double _dl2dx(Point2D point)
    => alphas[1, 1];

    private double _dl3dx(Point2D point)
    => alphas[2, 1];

    private double _dl1dy(Point2D point)
    => alphas[0, 2];

    private double _dl2dy(Point2D point)
    => alphas[1, 2];

    private double _dl3dy(Point2D point)
    => alphas[2, 2];

    private int FindElement(Point2D point)
    {
        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
            vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
            vertices[2] = grid.Nodes[grid.Elements[ielem][2]];

            double s01 = Math.Abs((vertices[1].X - vertices[0].X) * (point.Y - vertices[0].Y) -
                     (point.X - vertices[0].X) * (vertices[1].Y - vertices[0].Y));

            double s12 = Math.Abs((vertices[2].X - vertices[1].X) * (point.Y - vertices[1].Y) -
                     (point.X - vertices[1].X) * (vertices[2].Y - vertices[1].Y));

            double s20 = Math.Abs((vertices[0].X - vertices[2].X) * (point.Y - vertices[2].Y) -
                     (point.X - vertices[2].X) * (vertices[0].Y - vertices[2].Y));

            double dD = Math.Abs(DeterminantD());

            if (Math.Abs(dD - (s01 + s12 + s20)) < 1e-10)
                return ielem;
        }
        return -1;
    }
    public double ValueAtPoint(Point2D point)
    {
        double res = 0;

        int ielem = FindElement(point);

        if (ielem != -1)
        {
            CalcuclateAlphas();
            for (int i = 0; i < basis.Length; i++)
            {
                res += solution[grid.Elements[ielem][i]] * basis[i](point);
            }
        }
        return res;
    }

    public void PrintSolution()
    {
        for (int i = 0; i < solution.Length; i++)
        {
            Console.WriteLine($"x = {grid.Nodes[i].X:e3}, y = {grid.Nodes[i].Y:e3}, Az = {solution[i]}");
        }
    }

    public static Basis Mult(Basis fst, Basis scnd) 
        => (point) => fst(point) * scnd(point);
    public static Basis Sum(Basis fst, Basis scnd)
        => (point) => fst(point) + scnd(point);
}

