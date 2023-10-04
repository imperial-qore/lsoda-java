package odesolver;

import static odesolver.SupportingFunctions.*;
import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class SupportingFunctionsTest {

    @Test
    void calculate_vector_norm(){
        double[] v = {0, 3.0, 5.0, 7.0, 9.0};
        double[] w = {0, 0.4, 0.2, 0.05, 0.35};
        assertEquals(9*0.35, vmnorm(v.length-1, v, w));
    }

    @Test
    void calculate_matrix_norm(){
        double[][] a = {{0,0,0,0,0,0},
                {0,1,2,3,4,5},  // 10
                {0,-5,-4,-3,-2,-1}, // 20
                {0,-6,4,0,0,2}, // 30
                {0,3,4,6,-4,-10},   // 38
                {0,4,5,-100,-4,3}}; //44+5/6
        double[] w = {0, 0.1, 0.2, 0.3, 0.2, 0.1};
        assertTrue(Math.abs(fnorm(w.length-1, a, w)-(44+5.0/6))<1e-14);
    }

    @Test
    void calculate_error_tolerance(){
        double[] ycur = {0, 5, 7, 9, 11};
        double[] rtol = {0, 1e-7, 1e-8, 1e-4, 1e-6};
        double[] atol = {0, 1e-4, 1e-8, 1e-4, 1e-5};
        // scalar rtol and atol
        double[] res1 = {0,5e-7+1e-4,7e-7+1e-4,9e-7+1e-4,1.1e-6+1e-4};
        double[] ewt = SupportingFunctions.ewset(ycur, 1, rtol, atol,ycur.length-1);
        for (int i=1; i<ycur.length;i++)
            assertTrue(Math.abs(res1[i]-ewt[i])<1e-16);

        //  scalar rtol and array atol
        double[] res2 = {0,5e-7+1e-4,7.1e-7, 9e-7+1e-4, 1.11e-5};
        ewt = SupportingFunctions.ewset(ycur, 2, rtol, atol,ycur.length-1);
        for (int i=1; i<ycur.length;i++)
            assertTrue(Math.abs(res2[i]-ewt[i])<1e-16);

        // array rtol and scalar atol
        double[] res3 = {0, 5e-7+1e-4, 7e-8+1e-4, 1e-3, 1.11e-4};
        ewt = SupportingFunctions.ewset(ycur, 3, rtol, atol,ycur.length-1);
        for (int i=1; i<ycur.length;i++)
            assertTrue(Math.abs(res3[i]-ewt[i])<1e-16);

        // array rtol and atol
        double[] res4 = {0, 5e-7+1e-4, 8e-8, 1e-3, 2.1e-5};
        ewt = SupportingFunctions.ewset(ycur, 4, rtol, atol,ycur.length-1);
        for (int i=1; i<ycur.length;i++)
            assertTrue(Math.abs(res4[i]-ewt[i])<1e-16);
    }

    @Test
    void solve_linear_system_in_lsoda(){
        double[][] arr = {{0,0,0,0},{0,2,1,5},{0,4,4,-4},{0,1,3,1}};
        ReturningValues res = Utility.LUDecomposition(arr,3);
        double[] y = {0,5,0,6};
        double[] analyticalRes = {0, -1, 2, 1};
        // miter != 2
        assertArrayEquals(y, solsy(y, res.a, y.length-1, res.ipvt, 1));
        // miter == 2
        double[] x = solsy(y, res.a, y.length-1, res.ipvt, 2);
        for (int i=1; i<y.length; i++)
            assertTrue(Math.abs(analyticalRes[i]-x[i])<1e-15);
    }

    @Test
    void illegal_input_for_interpolate(){
        double ETA = 2.2204460492503131e-16;
        double[] dky = {0, 3, 3};
        double[][] yh = {{0,0,0}, {0, 3.0, 3.0},{0, 4.0, 4.0}};
        // illegal order
        double[] res1 = intdy(100.0, 4, dky, 2, 100.5, 0.2,
                ETA, 0.2, yh, 2);
        assertEquals(-1, res1[0]);

        // illegal interval
        double[] res2 = intdy(100.2, 1, dky, 2, 100.5, 0.2,
                ETA, 0.2, yh, 2);
        assertEquals(-2, res2[0]);
    }

    @Test
    void interpolate_and_calculate_y(){
        double ETA = 2.2204460492503131e-16;
        // using the example:
        // y = x*e^x + 3
        // y' = f = (1+x)*e^x
        // y'' = (2+x)*e^x
        // y''' = (3+x)*e^x
        // compute yh when tn=1.1, h=0.1
        double h=0.1;
        double[] y = {0, 1.1*Math.pow(Math.exp(1),1.1)+3};
        double[][] yh = {{0,0},
                {0, 1.1*Math.pow(Math.exp(1),1.1)+3},
                {0,h*2.1*Math.pow(Math.exp(1),1.1)},
                {0,h*h*3.1*Math.pow(Math.exp(1),1.1)/2},
                {0,Math.pow(h,3)*4.1*Math.pow(Math.exp(1),1.1)/6}};
        double[] res = intdy(1.05, 0, y, 3, 1.1, h, ETA,
                h, yh, 1);
        double truncErr = Math.pow(h,4)*5.1*Math.pow(Math.exp(1),1.1)/24;
        double exactSol = 1.05*Math.pow(Math.exp(1),1.05)+3;
        assertTrue(Math.abs(res[1]-exactSol) <= truncErr + Math.pow(h,5));
        System.out.printf("The truncation error is: %f when h=%f and order=%d\n",
                exactSol-res[1],h,3);
    }

    @Test
    void interpolate_and_calculate_first_order_derivative(){
        double ETA = 2.2204460492503131e-16;
        // using the example:
        // y' = (1+x)*e^x
        // y'' = (2+x)*e^x
        // y''' = (3+x)*e^x
        // compute yh when tn=1.1, h=0.1
        double h=0.1;
        double[] y = {0, 3.1*Math.pow(Math.exp(1),1.1)};
        double[][] yh = {{0,0},
                {0, 1.1*Math.pow(Math.exp(1),1.1)+3},
                {0,h*2.1*Math.pow(Math.exp(1),1.1)},
                {0,h*h*3.1*Math.pow(Math.exp(1),1.1)/2},
                {0,Math.pow(h,3)*4.1*Math.pow(Math.exp(1),1.1)/6},
                {0,Math.pow(h,4)*5.1*Math.pow(Math.exp(1),1.1)/24}};

        double[] res = intdy(1.05, 1, y, 4, 1.1, h, ETA,
                h, yh, 1);
        double truncErr = Math.pow(h,5)*6.1*Math.pow(Math.exp(1),1.1)/120;
        double exactSol = 2.05*Math.pow(Math.exp(1),1.05);
        assertTrue(Math.abs(res[1]-exactSol) <= truncErr+Math.pow(h,5));
    }
}