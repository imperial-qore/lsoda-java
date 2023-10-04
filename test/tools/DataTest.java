package tools;

import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

class DataTest {

    @Test
    void test_write_file() throws IOException {
        String filePath = "data";
        String fileName = "test.csv";
        Data data = new Data(filePath, fileName,2);
        double t = 1.0;
        double[] y = {0,1.0,2.0};
        String str = t+","+y[1]+","+y[2];
        data.write(t, y);

        // verify the written content
        BufferedReader bufferedReader = new BufferedReader(new FileReader(filePath+"/"+fileName));
        String header = bufferedReader.readLine();
        assertEquals(header, "t,y_1,y_2");
        String line = bufferedReader.readLine();
        assertEquals(str, line);
        data.closeWriter();
    }
}