package odesolver;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
public class ITaskAndContinuationTest {

    LSODA lsoda = new LSODA(0,0,1.0e-12, 1.0e-12,12, 5);
    FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
        @Override
        public int getDimension() {
            return 1;
        }

        @Override
        public void computeDerivatives(double v, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
            // y = e^{-15x}
            // f = -15*y
            ydot[0] = -15*y[0];
        }
    };
    private int neq = 1;
    private double[] y = {1};

    private double t = 0;
    private double tout = 1;
    private int itol = 1;
    private double[] atol = {0, 1.0e-12};
    private double[] rtol = {0, 1.0e-12};
    private int istate = 1;
    private int iopt=0;
    private double tcrit, h0, hmax, hmin=0;
    private int mxstep, mxhnil, mxordn, mxords=0;
    private int ixpr=1;

    @BeforeEach
    void initialize(){
        lsoda.ode = ode;
    }

    @Test
    void test_continuation_with_itask1(){
        int itask = 1;
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        double[] res = {0, Math.exp(-15*(tout))};
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
        // Continue
        istate = 2;
        lsoda.lsoda(neq,y,t,2,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        res[1] = Math.exp(-15*2);
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
        istate = 3;
        // Continue and change the params but tn has pass tout
        lsoda.lsoda(neq,y,t,2-1e-10,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, 3, 10, tcrit, h0, hmax, hmin);
        res[1] = Math.exp(-15*(2-1e-10));
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
    }

    @Test
    void test_advance_one_step_and_return(){
        // itask == 2
        int itask = 2;
        h0=1.3216372009101796E-7;
        mxstep= mxhnil= mxordn= mxords= 0;
        tcrit= hmin= 0.0;
        hmax = 1;
        ixpr=iopt=1;
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        double[] res = {0, Math.exp(-15*(t+h0))};
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
        // Continue and change params
        lsoda.lsoda(neq,y,t,h0-1e-10,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, 3, 10, tcrit, h0, hmax, hmin);
        res[1] = Math.exp(-15*(h0-1e-10));
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-5);
    }

    @Test
    void test_stop_at_the_first_internal_mesh_point(){
        // itask == 3
        int itask = 3;
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        double[] res = {0, Math.exp(-15*(lsoda.tn))};
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
        // Continue and change parameters
        istate = 3;
        lsoda.lsoda(neq,y,t,3,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, 3, 5, tcrit, h0, hmax, hmin);
        res[1] = Math.exp(-15*3);
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
    }

    @Test
    void test_compute_at_tout_and_not_overshoot_tcrit() throws IOException {
        int itask = 4;
        tcrit = 1+1e-5;
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        // 1. not overshoot t_crit
        if (lsoda.write){
            BufferedReader bufferedReader = new BufferedReader(new FileReader("data/lsoda.csv"));
            String line = bufferedReader.readLine();
            boolean flag = true;
            while((line=bufferedReader.readLine())!=null){
                double t = Double.parseDouble(line.split(",")[0]);
                if (tcrit<t){
                    flag = false;
                }
            }
            bufferedReader.close();
            assertTrue(flag);
        }
        // 2. test the result
        double[] res = {0, Math.exp(-15*(tout))};
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
        // Continue
        istate = 2;
        lsoda.lsoda(neq,y,t,tout+1e-8,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        res[1] = Math.exp(-15*(tout+1e-8));
        assertTrue(Math.abs(lsoda.y[1]-res[1])<1e-6);
    }

    @Test
    void test_advance_one_step_but_not_pass_tcrit(){
        int itask = 5;
        tcrit = 1+1e-3;
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        double[] res = {0, Math.exp(-15*(lsoda.tn))};
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
        // Continue
        istate = 3;
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        res[1] = Math.exp(-15*(lsoda.tn));
        assertTrue(Math.abs(res[1]-lsoda.y[1])<1e-8);
        assertTrue(lsoda.tn<=tcrit);

    }

}