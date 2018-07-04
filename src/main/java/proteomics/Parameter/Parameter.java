package proteomics.Parameter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.ECL2;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Parameter {

    private static final Logger logger = LoggerFactory.getLogger(Parameter.class);
    private static final Pattern commentLinePattern = Pattern.compile("^#.*");
    private static final Pattern linePattern = Pattern.compile("([^#]+)=([^#]+)#*.*");
    private static final Pattern enzymePattern = Pattern.compile("(.+)\\s+([01])\\s+([A-Z]+)\\s+([A-Z\\-]+)");
    private static final Pattern debugScanPattern = Pattern.compile("^[0-9]+$");

    private Map<String, String> parameterMap = new LinkedHashMap<>(); // keep the parameter order as it is in the parameter.def

    public Parameter(String parameterFile) throws Exception {
        BufferedReader parameterReader = new BufferedReader(new FileReader(parameterFile));
        String line = parameterReader.readLine().trim();
        if (!line.contentEquals("# " + ECL2.version)) {
            throw new Exception(String.format(Locale.US, "The parameter file version (%s) is not compatible with current ECL2 version (%s).", line.substring(2), ECL2.version));
        }
        while ((line = parameterReader.readLine()) != null) {
            line = line.trim();
            Matcher commentLineMatcher = commentLinePattern.matcher(line);
            if (!commentLineMatcher.matches()) {
                // This is not a comment line
                Matcher lineMatcher = linePattern.matcher(line);
                if (lineMatcher.matches()) {
                    String parameterName = lineMatcher.group(1).trim();
                    String parameterValue = lineMatcher.group(2).trim();
                    parameterMap.put(parameterName, parameterValue);
                } else {
                    Matcher debugScanMatcher = debugScanPattern.matcher(line);
                    if (debugScanMatcher.matches()) {
                        parameterMap.put(debugScanMatcher.group(), "");
                    } else {
                        Matcher enzymeMatcher = enzymePattern.matcher(line);
                        if (enzymeMatcher.matches()) {
                            if (parameterMap.containsKey("enzyme_name_1")) {
                                parameterMap.put("enzyme_name_2", enzymeMatcher.group(1).trim());
                                parameterMap.put("is_from_C_term_2", enzymeMatcher.group(2).trim());
                                parameterMap.put("cleavage_site_2", enzymeMatcher.group(3).trim());
                                parameterMap.put("protection_site_2", enzymeMatcher.group(4).trim());
                            } else {
                                parameterMap.put("enzyme_name_1", enzymeMatcher.group(1).trim());
                                parameterMap.put("is_from_C_term_1", enzymeMatcher.group(2).trim());
                                parameterMap.put("cleavage_site_1", enzymeMatcher.group(3).trim());
                                parameterMap.put("protection_site_1", enzymeMatcher.group(4).trim());
                            }
                        }
                    }
                }
            }
        }
        parameterReader.close();
    }


    public Map<String, String> returnParameterMap() {
        return parameterMap;
    }
}