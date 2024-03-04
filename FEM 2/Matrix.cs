namespace UMFCourseProject;

public class Matrix
{
   private double[,] mat;
   public int Size { get; }

   public double this[int i, int j]
   {
      get => mat[i, j];
      set => mat[i, j] = value;
   }

   public Matrix(int size)
   {
      mat = new double[size, size];
      Size = size;
   }

   public void Clear()
       => Array.Clear(mat, 0, mat.Length);

   public void Copy(Matrix destination)
   {
      for (int i = 0; i < destination.Size; i++)
      {
         for (int j = 0; j < destination.Size; j++)
         {
            destination[i, j] = mat[i, j];
         }
      }
   }

   public static Matrix operator +(Matrix fstMatrix, Matrix sndMatrix)
   {
      Matrix resultMatrix = new(fstMatrix.Size);

      for (int i = 0; i < resultMatrix.Size; i++)
      {
         for (int j = 0; j < resultMatrix.Size; j++)
         {
            resultMatrix[i, j] = fstMatrix[i, j] + sndMatrix[i, j];
         }
      }

      return resultMatrix;
   }
   
   public static Matrix operator *(double coef, Matrix Mat)
   {
      Matrix resultMatrix = new(Mat.Size);

      for (int i = 0; i < resultMatrix.Size; i++)
      {
         for (int j = 0; j < resultMatrix.Size; j++)
         {
            resultMatrix[i, j] = coef * Mat[i, j];
         }
      }

      return resultMatrix;
   }
}