namespace FEM2;
public class Vector
{
   private double[] vec;
   public int Length { get; init; }

   public Vector(int dim)
   {
      vec = new double[dim];
      Length = dim;
   }

   public double this[int index]
   {
      get => vec[index];
      set => vec[index] = value;
   }

   public void Fill(double value)
   {
      for (int i = 0; i < Length; i++)
      {
         vec[i] = value;
      }
   }

   public static void Copy(Vector source, Vector destination)
   {
      for (int i = 0; i < source.Length; i++)
         destination[i] = source[i];
   }

   public void Norming()
   {
      double norm = Norm();

      for (int i = 0; i < Length; i++)
         vec[i] /= norm;
   }

   public double Norm()
   {
      double result = 0;

      for (int i = 0; i < Length; i++)
         result += vec[i] * vec[i];

      return Math.Sqrt(result);
   }

   public static Vector operator *(Matrix matrix, Vector vector)
   {
      Vector result = new(vector.vec.Length);

      for (int i = 0; i < vector.Length; i++)
         for (int j = 0; j < vector.Length; j++)
            result.vec[i] += matrix[i, j] * vector.vec[j];

      return result;
   }

   public static Vector operator -(Vector fstVector, Vector sndVector)
   {
      Vector result = new(fstVector.Length);

      for (int i = 0; i < fstVector.Length; i++)
         result[i] = fstVector.vec[i] - sndVector.vec[i];

      return result;
   }

   public static Vector operator +(Vector fstVector, Vector sndVector)
   {
      Vector result = new(fstVector.Length);

      for (int i = 0; i < fstVector.Length; i++)
         result[i] = fstVector.vec[i] + sndVector.vec[i];

      return result;
   }

   public static Vector operator *(double coef, Vector vector)
   {
      Vector result = new(vector.Length);

      for (int i = 0; i < vector.Length; i++)
         result[i] = vector.vec[i] * coef;

      return result;
   }

   public static double operator *(Vector fstVector, Vector sndVector)
   {
      double result = 0;

      for (int i = 0; i < fstVector.Length; i++)
         result += fstVector.vec[i] * sndVector.vec[i];

      return result;
   }
}