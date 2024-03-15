namespace FEM2;

public class Solver
{
    public delegate void SolveFunc();

    public SolveFunc Solve;
    private SparseMatrix matrix = default!;
    private Matrix matrixx = default!;
    private Vector vector = default!;
    public Vector solution = default!;
    private double eps;
    private int maxIter;

    public Solver(int maxIter, double eps)
    {
        this.eps = eps;
        this.maxIter = maxIter;
    }

    public void SetSLAE(Vector vector, SparseMatrix matrix)
    {
        this.vector = vector;
        this.matrix = matrix;

        Solve = CGM;
    }

    public void SetParametres(int maxIter, double eps)
    {
        this.eps = eps;
        this.maxIter = maxIter;
    }


    public void SetSLAE(Vector vector, Matrix matrix)
    {
        this.vector = vector;
        matrixx = matrix;

        Solve = CGMMatrix;
    }

    public void CGM()
    {
        try
        {
            ArgumentNullException.ThrowIfNull(matrix, $"{nameof(matrix)} cant be NULL.");
            ArgumentNullException.ThrowIfNull(vector, $"{nameof(vector)} cant be NULL.");
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }

        double vectorNorm = vector.Norm();

        solution = new(vector.Length);
        Vector z = new(vector.Length);

        Vector r = vector - matrix * solution;
        Vector.Copy(r, z);

        int iter;

        for (iter = 0; iter < maxIter && r.Norm() / vectorNorm >= eps; iter++)
        {
            var tmp = matrix * z;
            var alpha = r * r / (tmp * z);
            solution += alpha * z;
            var squareNorm = r * r;
            r -= alpha * tmp;
            var beta = r * r / squareNorm;
            z = r + beta * z;
        }

        Console.WriteLine($"Last iteration - {iter}\n" +
           $"Residual norm - {r.Norm() / vectorNorm}");
    }

    public void CGMMatrix()
    {
        try
        {
            ArgumentNullException.ThrowIfNull(matrixx, $"{nameof(matrixx)} cant be NULL.");
            ArgumentNullException.ThrowIfNull(vector, $"{nameof(vector)} cant be NULL.");
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }

        double vectorNorm = vector.Norm();

        solution = new(vector.Length);
        Vector z = new(vector.Length);

        Vector r = vector - matrixx * solution;
        Vector.Copy(r, z);

        int iter;

        for (iter = 0; iter < maxIter && r.Norm() / vectorNorm >= eps; iter++)
        {
            var tmp = matrixx * z;
            var alpha = r * r / (tmp * z);
            solution += alpha * z;
            var squareNorm = r * r;
            r -= alpha * tmp;
            var beta = r * r / squareNorm;
            z = r + beta * z;
        }

        Console.WriteLine($"Last iteration - {iter}\n" +
           $"Residual norm - {r.Norm() / vectorNorm}");
    }

    public void PrintSolution()
    {
        for (int i = 0; i < solution.Length; i++)
        {
            Console.WriteLine(solution[i]);
        }
    }
}
