package odesolver;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

class UtilityTest {

    @Test
    void compute_dot_product() {
        double[] dx = {0, 3.3, 5.5, 7.7, 9.9, 11.11};
        double[] dy = {0, 2.2, 4.4, 6.6, 8.8, 12.12};
        // n<=0
        assertEquals(0, Utility.Dot(-1, dx, 1, dy, 1));
        // the increments are not equal
        assertEquals(dx[1]*dy[1]+dx[2]*dy[3],Utility.Dot(2, dx,1,dy,2));
        // the increments are negative
        assertEquals(dx[5]*dy[5]+dx[3]*dy[3]+dx[1]*dy[1],Utility.Dot(3, dx, -2, dy, -2));
        // the increments are unit
        double sum=0;
        for (int i=1; i<dx.length;i++){
            sum += dx[i]*dy[i];
        }
        assertEquals(sum,Utility.Dot(5, dx, 1, dy, 1));
        // the increments are euqal but not unit
        assertEquals(dx[1]*dy[1]+dx[3]*dy[3]+dx[5]*dy[5],Utility.Dot(3, dx, 2, dy, 2));
    }

    @Test
    void find_max_magnitude() {
        double[] dx = {0, 7.7, 9.9, 11.11, 3.3, 5.5};
        // n <= 0
        assertEquals(0, Utility.findMaxMagnitude(-2, dx, 1));
        // the increment is negative
        assertEquals(1, Utility.findMaxMagnitude(5, dx, -2));
        // increments are not unit
        assertEquals(2, Utility.findMaxMagnitude(3, dx, 2));
        // increments are unit
        assertEquals(3, Utility.findMaxMagnitude(3, dx, 1));
    }

    @Test
    void calculate_AX() {
        double[] dx = {0, 3.3, 5.5, 7.7, 9.9, 11.11, 13.13, 15.15, 17.17, 19.19};
        double a = 1.99;
        // illegal input
        assertEquals(dx, Utility.calAX(5, a, dx, 5, 1));
        // increments are unit
        double[] dy = new double[dx.length];
        for (int i=1; i<dx.length; i++)
            dy[i] = a*dx[i];
        assertArrayEquals(dy, Utility.calAX(dx.length-1, a, dx, 1, 1));
        // increments are not unit
        for (int i=1; i<dx.length; i=i+2)
            dy[i] *= a;
        assertArrayEquals(dy, Utility.calAX(dx.length/2, a, dx, 2, 1));
    }

    @Test
    void calculate_AX_plus_Y() {
        double[] dx = {0, 3.3, 5.5, 7.7, 9.9, 11.11};
        double[] dy = {0, 2.2, 4.4, 6.6, 8.8, 10.1};
        double a = 2.3;
        // illegal input
        assertArrayEquals(dy, Utility.calAXPlusY(4, a, dx, 4, dy, 5, 3, 5));
        // a=0
        assertArrayEquals(dy, Utility.calAXPlusY(dx.length-1, 0, dx, 1, dy, 1, 1, 1));
        // nonpositive increments
        double[] dz = new double[dx.length];
        for (int i=1; i<dx.length;i++)
            if (i%2==0)
                dz[i] = dy[i];
            else
                dz[i] = a*dx[i]+dy[i];
        assertArrayEquals(dz, Utility.calAXPlusY(dx.length/2, a, dx, -2, dy, -2, 1, 1));
        // nonequal increments
        for (int i=1; i<3;i++)
            dz[i]=a*dx[(i-1)*2+1]+dy[i];
        assertArrayEquals(dz, Utility.calAXPlusY(2, a, dx, 2, dy, 1, 1, 1));
        // increments are equal but not unit
        for (int i=1; i<dx.length;i=i+2)
            dz[i]=a*dx[i]+dy[i];
        assertArrayEquals(dz, Utility.calAXPlusY(dx.length/2, a, dx, 2, dy, 2, 1, 1));
        // increments are unit
        double[] dx2 = {0, 3.3, 5.5, 7.7, 9.9, 11.11, 13.13, 15.15, 17.17, 19.19};
        double[] dy2 = {0, 2.2, 4.4, 6.6, 8.8, 10.1, 12.12, 14.14, 16.16, 18.18};
        double[] dz2 = new double[dx2.length];
        dz2[1]=dy2[1];
        for (int i=2; i<dx2.length;i++)
            dz2[i]=a*dx2[i]+dy2[i];
        assertArrayEquals(dz2, Utility.calAXPlusY(dx2.length-2, a, dx2, 1, dy2, 1, 2, 2));
    }

    @Test
    void perform_LU_decomposition_with_singular_matrix() {
        // singularity
        double[][] arr1 = {{0,0,0,0},{0,2,1,5},{0,4,4,-4},{0,2,1,5}};
        ReturningValues res1 = Utility.LUDecomposition(arr1, 3);
        assertTrue(res1.info!=0);

        double[][] arr2 = {{0,0,0,0},{0,2,4,3},{0,0,0,0},{0,5,3,7},{0,2,0,1}};
        ReturningValues res2 = Utility.LUDecomposition(arr2, 3);
        assertTrue(res2.info!=0);
    }

    @Test
    void solve_linear_system() {
        double[][] arr = {{0,0,0,0},{0,1,-1,3,7},{0,10,-1,0,-2},{0,100,2,2,4},{0,5,99,2,9}};
        ReturningValues res = Utility.LUDecomposition(arr,4);
        double[] b = {0,727.4,-175,740.4,1471.2};
        b = Utility.solveLinearSys(res.a, 4, res.ipvt, b); //[3.1,5.4,9.2,100.3]
        double[] analyticalRes = {0, 3.1, 5.4, 9.2, 100.3};

        for (int i=1; i<b.length; i++)
            assertTrue(Math.abs(b[i]-analyticalRes[i])<1e-14);
    }
}