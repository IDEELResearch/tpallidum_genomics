
import java.io.*;
import java.util.*;

public class GetTPseq_fix_replaceDot {

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		BufferedReader br = null;
		BufferedWriter bw = null;

		String s = "";
		
		String sampleFile="/pine/scr/h/e/hennelly/unc_wgs-05/ss14_pipeline/consensus_generate/samples.txt";
		br=new BufferedReader(new FileReader(new File(sampleFile)));
		HashMap<String,String> samplemap=new HashMap<String,String>();
		ArrayList<String> sampleList = new ArrayList<String>();
		//	SRRXXX sample name (## see samples.txt for details) This step is used to ID convert.
		while((s=br.readLine())!=null){
			String tmp[]=s.split("\t");
			samplemap.put(tmp[0], tmp[1]);
			sampleList.add(tmp[0]);
		}
		br.close();
		System.out.println("sampleList.size "+sampleList.size());
		
		for (int sampleIndex = 0; sampleIndex < sampleList.size(); sampleIndex++) {
			String sample = sampleList.get(sampleIndex);
//SS14
			String vcfFilePath = "/pine/scr/h/e/hennelly/unc_wgs-05/ss14_pipeline/consensus_generate/SNPs_only.recode.vcf";
			String depthFilePath = "/pine/scr/h/e/hennelly/unc_wgs-05/ss14_pipeline/depth/" + sample + ".depth.txt";
				
			String refFilePath = "/proj/ideel/resources/genomes/Tpallidum/SS14_CP004011.1.fasta";
//

			
			
			br = new BufferedReader(new FileReader(new File(refFilePath)));
			HashMap<String, String> refseq = new HashMap<String, String>();
			int count = 1;
			while ((s = br.readLine()) != null) {
				if (!s.startsWith(">")) {
					char tmp[] = s.toCharArray();
					for (char c : tmp) {
						refseq.put(count + "", c + "");
						count++;
					}
				}
			}
			br.close();
			System.out.println("refseq.size  " + refseq.size());
			//depth ‘ if (Integer.parseInt(tmp[3]) > 0)’ You could change 0 to 2 for getting  3 X coverage
			br = new BufferedReader(new FileReader(new File(depthFilePath)));
			HashSet<String> coverSet = new HashSet<String>();
			while ((s = br.readLine()) != null) {
				String tmp[] = s.split("\t");
				if (Integer.parseInt(tmp[2]) > 0) {
					coverSet.add(tmp[1]);
				}
			}
			br.close();
			System.out.println("coverSet.size  " + coverSet.size());
			
		

			br = new BufferedReader(new FileReader(new File(vcfFilePath)));
			HashMap<String, HashMap<String, String>> vcfmap = new HashMap<String, HashMap<String, String>>();

			while ((s = br.readLine()) != null) {
				if (!s.startsWith("#")) {
					String tmp[] = s.split("\t");
					String ref = tmp[3];
					String alts[] = tmp[4].split(",");
					List<String> altList = new ArrayList<String>();
					altList.add(ref);
					for (int i = 0; i < alts.length; i++) {
						altList.add(alts[i]);
					}
					
					for(int i=9;i<sampleList.size()+9;i++){
						String sampleName=sampleList.get(i-9);
						if (vcfmap.containsKey(sampleName)) {
							HashMap<String, String> tmpmap = vcfmap.get(sampleName);
							String index = tmp[i].split(":")[0];
							if (index.equals(".")) {
								index = 0 + "";
							}
							String alt = altList.get(Integer.parseInt(index));
							if (alt.equals("*")) {
								alt = "N";
							}
							if(tmp[i].split(":")[0].equals(".")){
								alt = "N";
								index=1000+"";
							}
							if (!index.equals("0")) {
								tmpmap.put(tmp[1], altList.get(0).length() + "\t"
										+ alt);
								vcfmap.put(sampleName, tmpmap);
							}

						} else {
							HashMap<String, String> tmpmap = new HashMap<String, String>();
							String index = tmp[i].split(":")[0];
							if (index.equals(".")) {
								index = 0 + "";
							}
							String alt = altList.get(Integer.parseInt(index));
							if (alt.equals("*")) {
								alt = "N";
							}
							if(tmp[i].split(":")[0].equals(".")){
								alt = "N";
								index=1000+"";
							}
							if (!index.equals("0")) {
								tmpmap.put(tmp[1], altList.get(0).length() + "\t"
										+ alt);
								vcfmap.put(sampleName, tmpmap);
							}

						}
					}
				
				}
			}
			br.close();
			System.out.println("vcfmap.size  " + vcfmap.size());

			// substitution
			HashMap<String, String> resseq = new HashMap<String, String>();
			for (int i = 1; i < refseq.size() + 1; i++) {
				String base = refseq.get(i + "");
				if (!coverSet.contains(i + "")) {
					base = "N";
				}
				resseq.put(i + "", base);
			}
			System.out.println("resseq.size " + resseq.size());
			
			//TP0433 CP004010.2:461695-462464
			//TP0470 CP004010.2:498736-499845
			//TP0897 CP004010.2:975921-977441
			//TP0433 CP004011.1:461728-462497
			//TP0470 CP004011.1:499709-498879
			//TP0897 CP004011.1:977124-976691
			
		//Nichols
			/*
			for(int i=461695;i<462464+1;i++){
				resseq.put(i + "", "N");
			}
			for(int i=498736;i<499845+1;i++){
				resseq.put(i + "", "N");
			}
			for(int i=975921;i<977441+1;i++){
				resseq.put(i + "", "N");
			}
			*/
		
//			SS14
			/*
			for(int i=461728;i<462497+1;i++){
				resseq.put(i + "", "N");
			}
			for(int i=499709;i<498879+1;i++){
				resseq.put(i + "", "N");
			}
			for(int i=977124;i<976691+1;i++){
				resseq.put(i + "", "N");
			}
			*/
			
			
			// vcf # Set up the output PATH here.
			BufferedWriter tmpbw = new BufferedWriter(new FileWriter(new File("/pine/scr/h/e/hennelly/unc_wgs-05/ss14_pipeline"+"/unc_wgs-05_ss14.haplocaller.joint.snp1.fasta"),true));
			tmpbw.append(">"+(samplemap.get(sample)));tmpbw.newLine();tmpbw.flush();
			
				HashMap<String, String> altmap = vcfmap.get(sample);
				for (int i = 1; i <resseq.size() + 1;) {
					String base = resseq.get(i + "");
					int length = 1;
					if (!base.equals("N") & altmap.containsKey(i + "")) {
						base = altmap.get(i + "").split("\t")[1];
						length = Integer.parseInt(altmap.get(i + "")
								.split("\t")[0]);
					}
					tmpbw.append(base);
					i = i + length;
				}
				tmpbw.newLine();tmpbw.flush();
				tmpbw.close();
	
		/*
		  BufferedWriter tmpbw = new BufferedWriter(new FileWriter(new File(
		  "H:/proj/ideel/jonbparrlab/users/fnindo/test7/variants/depth/seq/"
		 +sample+".fasta")));
		 tmpbw.append(">"+sample);tmpbw.flush();tmpbw.newLine();
		 HashMap<String,String> altmap=vcfmap.get(sample); 
		 for (int i = 1; i < resseq.size()
		  + 1; ) { String base = resseq.get(i+""); int length=1;
		  if(!base.equals("N") & altmap.containsKey(i+"")){ base=altmap.get(i+"").split("\t")[1];
		  length=Integer.parseInt(altmap.get(i+"").split("\t")[0]); }
		  if(resseq.get(i+"").equals("N")){ base="N"; } tmpbw.append(base);
		  i=i+length; }
		  
		  
		  tmpbw.flush();tmpbw.newLine(); tmpbw.close();
		 */
		}
	}

	public static String getRev(String seq) {

		char tmpchars[] = seq.toCharArray();
		HashMap<String, String> charmap = new HashMap<String, String>();
		charmap.put("A", "T");
		charmap.put("C", "G");
		charmap.put("T", "A");
		charmap.put("G", "C");
		charmap.put("N", "N");
		String res = "";
		for (int i = 0; i < tmpchars.length; i++) {
			res = res + charmap.get(tmpchars[i] + "");
		}
		return res;

	}
}
