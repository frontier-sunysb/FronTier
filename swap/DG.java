
import org.jsoup.Connection.Response;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

import java.text.*;
import java.util.*;
import java.io.*;

public class DG {
	public static double getRealTimePrice(String stockname) 
	{
	
	     int i = 0;
	     String htmlbase = "http://finance.yahoo.com/q?s="; // base address from yahoo finance
   	     String price = null;
	     String htmlAdr = htmlbase + stockname;
	     try {
	     Document doc = Jsoup.connect(htmlAdr).timeout(3000).get();
	     Element elprice = null;
//	     Element eltime = null;	
	     String priceid = "yfs_l84_" + stockname.toLowerCase();
//	     String timeid = "yfs_t53_" + stock[i];
	     elprice = doc.getElementById(priceid);
	     price = elprice.text();
//	     lineInfo.append(eltime.text().split("\\s+")[0]);
	     }
	     catch(Exception e)
	     {
				
	     }
	     return Double.parseDouble(price);
	}
	
	
}

