package edu.scripps.dia;

import java.util.ArrayList;
import java.util.List;

import static rpark.statistics.Outlier.removeOutlier;

/**
 * Created by Titus Jung titusj@scripps.edu on 7/31/19.
 */
public class Scrath {

    public static void main(String [] args){

        List<Double> test = new ArrayList<>();
      /*  test.add(124.14);
        test.add(124.24);
        test.add(124.34);
        test.add(124.46);
        test.add(124.57);
        test.add(1.0);*/
        test.add(95.8499984741211);
        test.add(47.0999984741211);
        test.add(52.0900001525879);
        test.add(69.4300003051758);
        test.add(55.7700004577637);
        test.add(32.189998626709);
        test.add(68.7600021362305);
        test.add(85.5500030517578);
        test.add(74.9300003051758);
        test.add(75.1999969482422);
        test.add(46.8400001525879);
        test.add(49.560001373291);
        test.add(110.809997558594);
        test.add(77.2399978637695);
        test.add(81.8199996948242);
        test.add(49.3199996948242);
        test.add(47.9199981689453);
        test.add(48.2200012207031);
        test.add(50.7700004577637);
        test.add(46.4300003051758);
        test.add(45.6100006103516);
        test.add(47.8300018310547);
        test.add(60.1199989318848);
        test.add(120.870002746582);
        test.add(7.80000019073486);
        test.add(50.2700004577637);
        test.add(51.7599983215332);
        test.add(51.6399993896484);
        test.add(0.56999999284744);
        test.add(66.9800033569336);
        test.add(102.139999389648);
        test.add(44.0299987792969);
        test.add(45.6699981689453);
        test.add(48.3400001525879);
        test.add(48.9799995422363);
        test.add(52.5099983215332);
        test.add(45.3300018310547);
        test.add(0.75);
        test.add(47.3600006103516);
        test.add(40.5200004577637);
        test.add(46.4599990844727);
        test.add(44.9099998474121);



        List<Double> filteredList = removeOutlier(test,0.1);
        System.out.println(test.size()+"\t"+filteredList.size());
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;
        for(double d: filteredList)
        {
            max = max<d ? d : max;
            min = min>d ? d : min;
        }
        System.out.println("max: "+max + "\tmin: "+min);

    }
}
