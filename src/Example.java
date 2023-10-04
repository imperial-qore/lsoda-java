import odesolver.LSODA;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import java.util.Arrays;

public class Example{

    public static void main(String[] args){
        LSODA lsoda = new LSODA(0,0,1.0e-8,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 3;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = 1.0e4 * y[1] * y[2] - 0.04 * y[0] ;
                ydot[2] = 3.0e7 * y[1] * y[1];
                ydot[1] = -1.0 * (ydot[0] + ydot[2]);
            }
        };
        double[] y = {1.0,0.0,0.0};
        double t = 0;
        double tout = 4e5;
        double[] result = new double[3];
        lsoda.integrate(ode,t,y,tout,result);
        System.out.println(Arrays.toString(result));
    }
}
