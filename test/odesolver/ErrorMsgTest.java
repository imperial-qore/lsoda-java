package odesolver;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
public class ErrorMsgTest {
    LSODA lsoda = new LSODA(0,0,1.0e-12,1.0e-12,12, 5);
    FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
        @Override
        public int getDimension() {
            return 1;
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
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
    private int istate = 1, ixpr=1, iopt=0,itask=1,mxstep, mxhnil, mxordn, mxords;
    private double tcrit, h0, hmax, hmin;
    @BeforeEach
    void initialize(){
        lsoda.ode = ode;
    }

    // block a
    @Test
    void test_wrong_itask(){
        itask = 0;
        try {
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(0.0, Math.abs(t - lsoda.tn));
    }

    @Test
    void test_wrong_istate(){
        istate = 4;
        try{
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println("RuntimeException");
        }
        assertEquals(0.0, Math.abs(t - lsoda.tn));
        istate = 2;
        try {
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(0.0, Math.abs(t - lsoda.tn));
    }

    @Test
    void test_repeated_input(){
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);

        assertTrue(Math.abs(Math.exp(-15*tout)-lsoda.y[1])<1e-8);
        try{
            for (int i=0; i<5;i++)
                lsoda.lsoda(neq,y,0,0,itol,rtol,atol,itask,1,iopt,
                        ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }


    }

    // block b
    @Test
    void test_number_of_equations(){
        try {
            lsoda.lsoda(-1,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(0.0, Math.abs(t - lsoda.tn));
        iopt=1;
        // change the number of equations after a successful step
        itask = 2;
        h0=1e-8;
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        istate = 3;
        try{
            lsoda.lsoda(5,y,t,t,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(0.0, Math.abs(h0 - lsoda.tn));
    }

    @Test
    void test_itol_iopt_and_optional_params(){
        try{
            lsoda.lsoda(neq,y,t,tout,0,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        try{
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,3,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(0.0, Math.abs(t - lsoda.tn));
        // optional params
        iopt = 1;
        try {
            // ixpr
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    -1, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }

        try{
            // max step
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, -5, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        try{
            // mxhnil
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, -10, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }

        try{
            // order
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, -10, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        try{
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, -9, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        try {
            // step size
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, -0.1, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }

        try {
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, -0.1);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        try {
            // tolerance
            rtol[1]=-1e-4;
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        try{
            atol[1]=-1e-4;
            rtol[1]= 1e-12;
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
    }

    // block c
    @Test
    void test_tcrit_and_tout(){
        try {
            itask = 4;
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tout-1e-5, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
    }

    @Test
    void test_initial_step_size(){
        // adjust step size according to the h_max
        itask = 2;
        hmax=1e-20;
        iopt = 1;
        try {
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tout, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }

        assertEquals(lsoda.tn,hmax);
    }

    // block d
    @Test
    void test_interpolation_flag_in_lsoda(){
        // itask = 1
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        // continuation
        double tcur=lsoda.tn;
        istate = 3;
        try{
            lsoda.lsoda(neq,y,t,tout-0.15,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(lsoda.tn, tcur);
        // itask = 4
        itask = 4;
        try {
            lsoda.lsoda(neq,y,t,tout-0.15,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tout+0.1, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(lsoda.tn, tcur);
    }

    @Test
    void test_position_check_among_t(){
        lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        double tcur = lsoda.tn;
        itask = 3;
        istate = 3;
        try{
            lsoda.lsoda(neq,y,t,tout-0.15,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(lsoda.tn, tcur);

        itask = 4;
        //tn and tout
        try {
            lsoda.lsoda(neq,y,t,tout-0.15,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(lsoda.tn, tcur);

        itask = 5;
        try {
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
        assertEquals(lsoda.tn, tcur);
    }

    // block e and f
    @Test
    void test_max_steps_taken(){
        iopt = 1;
        mxstep = 10;
        try {
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
    }

    @Test
    void test_much_accuracy_warning_in_operation(){
        atol[1] = 1.0e-15;
        rtol[1] = 1.0e-14;
        try{
            lsoda.lsoda(neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
                    ixpr, mxstep, mxhnil, mxordn, mxords, tcrit, h0, hmax, hmin);
        }catch (RuntimeException e){
            System.out.println(e.getMessage());
        }
    }
}
