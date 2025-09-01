import java.util.*;
import org.json.JSONObject;
import org.json.JSONException;
import java.io.BufferedReader;
import java.io.InputStreamReader;

public class HashiraPlacements {

    // Convert a number from any base to base 10
    public static long convertToBase10(String value, int base) {
        if (base == 10 || base == 0) { // base 0 treated as decimal
            return Long.parseLong(value);
        }

        long result = 0;
        String digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        for (int i = 0; i < value.length(); i++) {
            char c = value.charAt(i);
            int digit = digits.indexOf(Character.toUpperCase(c));
            if (digit < 0 || digit >= base) {
                throw new IllegalArgumentException("Invalid digit '" + c + "' for base " + base);
            }
            result = result * base + digit;
        }

        return result;
    }

    // Gaussian elimination to solve Ax = b
    public static double[] gaussianElimination(double[][] A, double[] b) {
        int n = b.length;

        for (int p = 0; p < n; p++) {
            // Find pivot row and swap
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }

            // Swap rows in A
            double[] tempA = A[p];
            A[p] = A[max];
            A[max] = tempA;

            // Swap in b
            double tempB = b[p];
            b[p] = b[max];
            b[max] = tempB;

            // Check for singularity
            if (Math.abs(A[p][p]) <= 1e-12) {
                throw new RuntimeException("Matrix is singular or nearly singular");
            }

            // Eliminate below pivot
            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // Back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        return x;
    }

    public static void main(String[] args) {
        try {
            // Read JSON input
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            StringBuilder jsonInput = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                jsonInput.append(line);
            }

            JSONObject json = new JSONObject(jsonInput.toString());

            // Extract k
            JSONObject key = json.getJSONObject("key");
            int k = key.getInt("k"); // number of coefficients to find

            // Collect (x, f(x)) points
            List<Long> xValues = new ArrayList<>();
            List<Long> fValues = new ArrayList<>();

            for (String jsonKey : json.keySet()) {
                if (!jsonKey.equals("key")) {
                    try {
                        JSONObject root = json.getJSONObject(jsonKey);
                        int base = Integer.parseInt(root.getString("base"));
                        String value = root.getString("value");

                        long decimalValue = convertToBase10(value, base);
                        xValues.add(decimalValue);
                        fValues.add(0L); // root => f(x) = 0
                    } catch (JSONException e) {
                        System.err.println("Skipping invalid entry: " + jsonKey);
                    }
                }
            }

            int n = xValues.size();

            // Validation
            if (n < k) {
                throw new RuntimeException("Not enough points to solve for " + k + " coefficients");
            }

            // Build matrix system
            double[][] A = new double[n][k];
            double[] b = new double[n];

            for (int i = 0; i < n; i++) {
                long x = xValues.get(i);
                for (int j = 0; j < k; j++) {
                    A[i][j] = Math.pow(x, j);
                }
                b[i] = fValues.get(i);
            }

            // Solve for coefficients
            double[] coefficients = gaussianElimination(A, b);

            // Print result
            System.out.println("Polynomial Coefficients:");
            for (int i = 0; i < coefficients.length; i++) {
                System.out.printf("a%d = %.2f%n", i, coefficients[i]);
            }

        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }
    }
}
