package odesolver;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class StiffExamplesTest {
    @Test
    void solve_van_der_pol_oscillator(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                int mu = 1000;
                ydot[0] = y[1] ;
                ydot[1] = mu*(1-y[0]*y[0])*y[1]-y[0];
            }
        };

        double[] y = {2.0,0.0};
        double t = 0;
        double tout = 3000;
        long startTime = System.nanoTime();
        double[] result = new double[2];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("y1 = %16.16e, y2 = %16.16e\n", result[0], result[1]);
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
    }

    /*
    * Class A: Linear with real eigenvalues
    * */
    @Test
    void ExampleA1(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -0.5*y[0] ;
                ydot[1] = -y[1];
                ydot[2] = -100 * y[2];
                ydot[3] = -90 * y[3];
            }
        };
        double[] y = {0.0,0.0,0.0,0.0};
        double t = 0;
        double tout = 120;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("y1 = %16.16e, y2 = %16.16e, y3 = %16.16e, y4 = %16.16e\n", result[0], result[1],result[2], result[3]);
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
    }

    @Test
    void ExampleA2(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 9;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -1800*y[0] + 900*y[1];
                for (int i=1;i<8;i++)
                    ydot[i] = y[i-1] - 2*y[i] + y[i+1];
                ydot[8] = 1000 * y[7] - 2000*y[8] + 1000;
            }
        };
        double[] y = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        double t = 0;
        double tout = 120;
        long startTime = System.nanoTime();
        double[] result = new double[9];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        for (int i=0; i<9; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
    }

    @Test
    void ExampleA3(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -10000 * y[0] + 100*y[1] - 10*y[2] + y[3];
                ydot[1] = -1000*y[1] + 10*y[2] - 10*y[3];
                ydot[2] = -y[2] + 10*y[3];
                ydot[3] = -0.1 * y[3];
            }
        };
        double[] y = {1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("y1 = %16.16e, y2 = %16.16e, y3 = %16.16e, y4 = %16.16e\n", result[0], result[1],result[2], result[3]);
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
    }

    @Test
    void ExampleA4(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 10;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                for (int i=0; i<10; i++)
                    ydot[i] = -Math.pow(i+1,5)*y[i];
            }
        };
        double[] y = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 1;
        long startTime = System.nanoTime();
        double[] result = new double[10];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        for (int i=0; i<10; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
    }


    /*
    * Class B: Linear with non-real eigenvalues
    * */
    @Test
    void ExampleB1(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -y[0] + y[1];
                ydot[1] = -100*y[0] - y[1];
                ydot[2] = -100*y[2] + y[3];
                ydot[3] = -10000 * y[2] - 100 * y[3];
            }
        };
        double[] y = {1.0,0,1.0,0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleB2(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 6;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                int a = 3;
                ydot[0] = -y[0] + a*y[1];
                ydot[1] = -a*y[0]-10*y[1];
                ydot[2] = -4 * y[2];
                ydot[3] = -y[3];
                ydot[4] = -0.5 * y[4];
                ydot[5] = -0.1 * y[5];
            }
        };
        double[] y = {1.0,1.0,1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[6];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<6; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleB3(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 6;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                int a = 8;
                ydot[0] = -y[0] + a*y[1];
                ydot[1] = -a*y[0]-10*y[1];
                ydot[2] = -4 * y[2];
                ydot[3] = -y[3];
                ydot[4] = -0.5 * y[4];
                ydot[5] = -0.1 * y[5];
            }
        };
        double[] y = {1.0,1.0,1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[6];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<6; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }
    @Test
    void ExampleB4(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 6;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                int a = 25;
                ydot[0] = -y[0] + a*y[1];
                ydot[1] = -a*y[0]-10*y[1];
                ydot[2] = -4 * y[2];
                ydot[3] = -y[3];
                ydot[4] = -0.5 * y[4];
                ydot[5] = -0.1 * y[5];
            }
        };
        double[] y = {1.0,1.0,1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[6];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<6; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleB5(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 6;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                int a = 100;
                ydot[0] = -y[0] + a*y[1];
                ydot[1] = -a*y[0]-10*y[1];
                ydot[2] = -4 * y[2];
                ydot[3] = -y[3];
                ydot[4] = -0.5 * y[4];
                ydot[5] = -0.1 * y[5];
            }
        };
        double[] y = {1.0,1.0,1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[6];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<6; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    /*
     * Class C: Non-linear Coupling
     * */

    @Test
    void ExampleC1(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3];
                ydot[1] = -10 * y[1] + 10 * (y[2]*y[2] + y[3]*y[3]);
                ydot[2] = -40 * y[2] + 40 * y[3]*y[3];
                ydot[3] = -100 * y[3] + 2;
            }
        };
        double[] y = {1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[6];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleC2(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                double b = 0.1;
                ydot[0] = -y[0]+ 2;
                ydot[1] = -10*y[1] + b*y[0]*y[0];
                ydot[2] = -40*y[2] + 40*b*(y[0]*y[0]+y[1]*y[1]);
                ydot[3] = -100*y[3] + 10*b*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
            }
        };
        double[] y = {1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[6];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }
    @Test
    void ExampleC3(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                double b = 1;
                ydot[0] = -y[0]+ 2;
                ydot[1] = -10*y[1] + b*y[0]*y[0];
                ydot[2] = -40*y[2] + 40*b*(y[0]*y[0]+y[1]*y[1]);
                ydot[3] = -100*y[3] + 10*b*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
            }
        };
        double[] y = {1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i+1,result[i]);
    }
    @Test
    void ExampleC4(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                double b = 10;
                ydot[0] = -y[0]+ 2;
                ydot[1] = -10*y[1] + b*y[0]*y[0];
                ydot[2] = -40*y[2] + 40*b*(y[0]*y[0]+y[1]*y[1]);
                ydot[3] = -100*y[3] + 10*b*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
            }
        };
        double[] y = {1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i+1,result[i]);
    }
    @Test
    void ExampleC5(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                double b = 20;
                ydot[0] = -y[0]+ 2;
                ydot[1] = -10*y[1] + b*y[0]*y[0];
                ydot[2] = -40*y[2] + 40*b*(y[0]*y[0]+y[1]*y[1]);
                ydot[3] = -100*y[3] + 10*b*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
            }
        };
        double[] y = {1.0,1.0,1.0,1.0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i+1,result[i]);
    }

    /*
    Class D: Non-linear with real eigenvalues
     */
    @Test
    void ExampleD1(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 3;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = 0.2 * (y[1]-y[0]);
                ydot[1] = 10*y[0] - (60 - 0.125 * y[2]) * y[1]+ 0.125 * y[2];
                ydot[2] = 1;
            }
        };
        double[] y = {0,0,0};
        double t = 0;
        double tout = 400;
        long startTime = System.nanoTime();
        double[] result = new double[3];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<3; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleD2(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 3;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -0.04 * y[0] + 0.01 * y[1] * y[2];
                ydot[1] = 400 * y[0] - 100*y[1]*y[2] - 3000 * y[1] * y[1];
                ydot[2] = 30 * y[1] * y[1];
            }
        };
        double[] y = {1.0, 0, 0};
        double t = 0;
        double tout = 40;
        long startTime = System.nanoTime();
        double[] result = new double[3];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<3; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleD3(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = y[2] - 100 * y[0] * y[1];
                ydot[1] = y[2] + 2 * y[3] - 100*y[0]*y[1] - 20000 * y[1] * y[1];
                ydot[2] = -y[2] + 100 * y[0] * y[1];
                ydot[3] = -y[3] + 10000 * y[1] * y[1];
            }
        };
        double[] y = {1.0, 1.0, 0, 0};
        double t = 0;
        double tout = 20;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleD4(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 3;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -0.013 * y[0] - 1000 * y[0] * y[2];
                ydot[1] = -2500 * y[1] * y[2];
                ydot[2] = -0.013 * y[0] - 1000 * y[0] * y[2] - 2500 * y[1] * y[2];
            }
        };
        double[] y = {1.0, 1.0, 0};
        double t = 0;
        double tout = 50;
        long startTime = System.nanoTime();
        double[] result = new double[3];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<3; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleD5(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 2;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = 0.01 - (1 + (y[0] + 1000)*(y[0] + 1))*(0.01 + y[0] + y[1]);
                ydot[1] = 0.01 - (1 + y[1]*y[1])*(0.01 + y[0] + y[1]);
            }
        };
        double[] y = {0, 0};
        double t = 0;
        double tout = 100;
        long startTime = System.nanoTime();
        double[] result = new double[2];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<2; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleD6(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 3;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -y[0] + Math.pow(10,8) * y[2] * (1 - y[0]);
                ydot[1] = -10 * y[1] + 3 * Math.pow(10,7) * y[2] * (1-y[1]);
                ydot[2] = -ydot[0] - ydot[1];
            }
        };
        double[] y = {1.0, 0, 0};
        double t = 0;
        double tout = 1;
        long startTime = System.nanoTime();
        double[] result = new double[3];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<3; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    /*
    Class E: Non-linear with non-real eigenvalues
     */
    @Test
    void ExampleE1(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }
            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                int gamma = 100;
                ydot[0] = y[1];
                ydot[1] = y[2];
                ydot[2] = y[3];
                ydot[3] = (y[0] * y[0] - Math.sin(y[0]) - Math.pow(gamma,4)) * y[0] +
                        (y[1] * y[2]/(y[0]*y[0]+1) - 4 * Math.pow(gamma,3)) * y[1] +
                        (1 - 6 * gamma*gamma) * y[2] + (10 * Math.exp(-y[3]*y[3]) - 4*gamma)*y[3] + 1;
            }
        };
        double[] y = {0, 0, 0, 0};
        double t = 0;
        double tout = 1;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleE2(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                int mu = 5;
                ydot[0] = y[1] ;
                ydot[1] = mu*(1-y[0]*y[0])*y[1]-y[0];
            }
        };

        double[] y = {2.0,0.0};
        double t = 0;
        double tout = 1;
        long startTime = System.nanoTime();
        double[] result = new double[2];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
        System.out.printf("y1 = %16.16e, y2 = %16.16e\n", result[0], result[1]);

    }

    @Test
    void ExampleE3(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 3;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -(55 + y[2])*y[0] + 65 * y[1];
                ydot[1] = 0.0785 * (y[0] - y[1]);
                ydot[2] = 0.1 * y[0];
            }
        };

        double[] y = {1.0,1.0,0.0};
        double t = 0;
        double tout = 500;
        long startTime = System.nanoTime();
        double[] result = new double[3];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
        for (int i=0; i<3; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }

    @Test
    void ExampleE4(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                double[][] U = {{-0.5, 0.5, 0.5, 0.5},
                        {0.5, -0.5, 0.5, 0.5},
                        {0.5, 0.5, -0.5, 0.5},
                        {0.5, 0.5, 0.5, -0.5}};
                double[] Z = new double[4];
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        Z[i] += U[i][j] * y[j];

                double[] Gy = new double[4];
                double[] tmp = {(Z[0] * Z[0] - Z[1] * Z[1]) / 2, Z[0] * Z[1], Z[2] * Z[2], Z[3] * Z[3]};
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        Gy[i] += U[i][j] * tmp[j];

                double[][] coe = {{-10, -10, 0, 0},
                        {10, -10, 0, 0},
                        {0, 0, 1000, 0},
                        {0, 0, 0, 0.01}};
                double[][] UeU = multiply(transpose(U),coe,4);
                UeU = multiply(UeU,U,4);
                for (int i=0; i<4;i++){
                    for (int j=0; j<4;j++)
                        ydot[i] -= UeU[i][j]*y[j];
                    ydot[i] += Gy[i];
                }

            }
            public double[][] multiply(double[][] a,double[][] b,int n){
                double[][] result = new double[n][n];
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++)
                        for(int k=0;k<n;k++)
                            result[i][j] += a[i][k]*b[k][j];
                }
                return result;
            }

            public double[][] transpose(double[][] U){
                double[][] result = new double[U.length][U.length];
                for(int i=0; i<U.length; i++) {
                    for(int j=0; j<U[i].length; j++) {
                        result[j][i] = U[i][j];
                    }
                }
                return result;
            }
        };

        double[] y = {0,-2,-1,-1};
        double t = 0;
        double tout = 500;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d\n", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }
    @Test
    void ExampleE5(){
        LSODA lsoda = new LSODA(0,0,1.0e-6,1.0e-6,12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 4;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                ydot[0] = -7.89 * Math.pow(10,-10)*y[0] - 1.1 * Math.pow(10,7)*y[0]*y[2];
                ydot[1] = 7.89 * Math.pow(10, -10)*y[0] - 1.13 * Math.pow(10,9)*y[1]*y[2];
                ydot[2] = 7.89 * Math.pow(10, -10)*y[0] - 1.1 * Math.pow(10,7) * y[0] * y[2] +
                        1.13 * Math.pow(10,3) * y[3] -1.13 * Math.pow(10,9)*y[1]*y[2];
                ydot[3] = 1.1 * Math.pow(10,7) * y[0] * y[2] - 1.13 * Math.pow(10,3)*y[3];
            }
        };

        double[] y = {1.76 * Math.pow(10,-3),0,0,0};
        double t = 0;
        double tout = 1000;
        long startTime = System.nanoTime();
        double[] result = new double[4];
        lsoda.integrate(ode,t,y,tout,result);
        long endTime = System.nanoTime();
        System.out.printf("Running time: %.5f ms\n", (endTime-startTime)/1000000.0);
        System.out.printf("number of Jacobian evals:%d\nnumber of f evals:%d\n",
                lsoda.getJacobianEvaluations(),lsoda.getEvaluations());
        System.out.printf("Steps taken: %d", lsoda.getStepsTaken());
        for (int i=0; i<4; i++)
            System.out.printf("y%d = %16.16e\n", i,result[i]);
    }
}
