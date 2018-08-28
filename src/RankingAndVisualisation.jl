include("SubstructureRanking.jl")
include("StructureVisualisation.jl")
using CSV

function rankbycoevolution(outputprefix, alignmentfile, maxfile, paired::Array{Int,1}, mapping, dataset::Dataset, lambdas::Array{Float64,1})
    outputdir = joinpath(dirname(outputprefix), "ranking/")
    if !isdir(outputdir)
        mkpath(outputdir)
    end
    ret, finalls, texoutput = substructureranking(paired, lambdas, AbstractString[], mapping, nothing, false, false, false)

    try
        template = read(open(joinpath(@__DIR__,"ranking_template.tex"),"r"), String)
        output = replace(template, "#INSERT#" => texoutput)
        texfile = abspath(joinpath(outputdir,"ranking.tex"))
        println("Texfile: ", texfile)
        fout = open(texfile,"w")
        write(fout, output)
        close(fout)


        run(`pdflatex -aux-directory="$(abspath(joinpath(outputdir)))" -output-directory="$(abspath(joinpath(outputdir)))" $(texfile)`)
    catch e

    end

    #sequence, dummy = readctfile("/media/michael/Sandisk500GB/data/hiv1-SHAPE-MaP.ct")

    #ctfile = "D:\\Dropbox\\dev\\farce-julia\\src\\mcmc\\hiv1b\\hiv1-SHAPE-MaP.ct"
    #csvfile = "D:\\Dropbox\\Final\\thesis\\thesis\\manuscript\\datasets\\NL4-3 SHAPE reactivities.csv"

    ctfile = "D:\\Dropbox\\dev\\farce-julia\\src\\mcmc\\hiv1b\\hiv1-SHAPE-MaP.ct"
    csvfile = "D:\\Dropbox\\Final\\thesis\\thesis\\manuscript\\datasets\\NL4-3 SHAPE reactivities.csv"

    #shapereactivities,shapeprobs,mapping, revmapping = loadNL43ShapeReactivities2(dataset, alignmentfile, ctfile, csvfile)
    #backdata = shapereactivities
    mapping = Int[i for i=1:length(paired)]
    backdata = zeros(Float64, length(paired))

    rank = 1
    for sub in finalls
        startpos = sub[3]
        endpos = sub[4]
        len = endpos - startpos + 1
        substructure = paired[startpos:endpos]
        for i=1:length(substructure)
            if substructure[i] > 0
                substructure[i] = substructure[i] - startpos + 1
                if substructure[i] < 0 || substructure[i] > len
                    substructure[i] = 0
                end
            end
        end

        data = zeros(Float64, size(dataset.data,1), len, 4)
        for i=1:size(data,1)
            for j=1:len
                for k=1:4
                    data[i,j,k] = dataset.data[i, j+startpos-1,k]
                end
            end
        end

        drawsubstructure(startpos::Int, dataset, alignmentfile, joinpath(outputdir, string(startpos,"-",endpos,".rank", rank,"of",length(finalls),".svg")), data, dataset.gapfrequency[startpos:endpos], substructure, backdata[startpos:endpos], lambdas[startpos:endpos], zeros(Float64,length(paired))[startpos:endpos], mapping)
        rank += 1
    end
end

function drawsubstructure(startpos::Int, dataset, alignmentfile, outfile, data::Array{Float64,3}, gapfrequency::Array{Float64,1}, paired::Array{Int,1}, backvalues::Array{Float64,1}, linkvalues::Array{Float64,1}, linkvalues2::Array{Float64,1}, mapping, visible::Array{Int,1}=Int[])
  maxcoevolutionvalue = 10.0

  ylgnbu = reverse([[1.0,1.0,0.85],[0.93,0.98,0.71],[0.78,0.91,0.71],[0.50,0.80,0.73],[0.27,0.71,0.76],[0.15,0.56,0.74],[0.15,0.37,0.64],[0.16,0.22,0.56],[0.06,0.14,0.34]])

  shapecolors = [(0,0,0,1.0),(0,0,0,1.0),  (33,177,79,1.0),(33,177,79,1.0),  (250,200,36,1.0),(250,200,36,1.0), (244,23,28,1.0),(244,23,28,1.0),(244,23,28,1.0)]
  shapepositions = [0.0, 0.3, 0.3, 0.4, 0.4, 0.6, 0.6, 1.0]
  shapegradient = ColorGradient(shapecolors, shapepositions)

  #linkgradient =  viridisgradient(ylgnbu)
  linkgradient =  viridisgradient(_viridis_data)
  backgradient =  ColorGradient([(255,255,255,1.0),(255,255,255,1.0)], [0.0,1.0])
  backlabel = ""

  #drawlegend(string(outfile, ".coevolutionprobability.svg"), linkgradient, ["0.0", "2.5", "5.0", "7.5", "10.0"])
  #exit()
  #backgradient = shapegradient
  #backlabel = "SHAPE reactivity"

  #=
  drawlegend(string(outfile, ".substitutionrate.svg"), backgradient, ["low", "moderate", "high"])
  drawlegend(string(outfile, ".coevolutionprobability.svg"), linkgradient, ["0.00", "0.25", "0.50", "0.75", "1.00"])
  drawlegend(string(outfile, ".shape.svg"), shapegradient, ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])=#
  points, minx, miny, maxx, maxy = get2Dpoints(paired)
  numpoints = length(points)

  scale = 2.0
  radius = 7.5*scale
  xoffset = radius
  yoffset = radius
  width = (maxx-minx)*scale + xoffset*2.0 + 50.0
  height = (maxy-miny)*scale + yoffset*2.0 + 50.0
  #width = 3000
  #height = 1500

  for i=1:numpoints
    points[i][1] = (points[i][1] - minx + xoffset)*scale
    points[i][2] = (points[i][2] - miny + yoffset)*scale
  end

  svgout = open(outfile,"w")
  write(svgout,"""
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" onload="init(evt)" width="$width" height="$height"  style=\"font-family:sans-serif\">
  """
  )
  write(svgout, """
  <defs>
  <marker id="markerSquare" markerWidth="7" markerHeight="7" refX="4" refY="4"
          orient="auto">
      <rect x="1" y="1" width="5" height="5" style="stroke: none; fill:#000000;"/>
  </marker>

  <marker id="markerArrow" markerWidth="13" markerHeight="13" refX="2" refY="7"
          orient="auto">
      <path d="M2,2 L2,13 L8,7 L2,2" style="fill: #000000;" />
  </marker>
</defs>

  <style>
    .caption{
	font-size: 14px;
	font-family: Georgia, serif;
    }
    .tooltip{
	font-size: 12px;
	color: black;
    }
    .tooltip_bg{
	fill: white;
	stroke: black;
	stroke-width: 1;
	opacity: 1.0;
    }
  </style>
  """)
  write(svgout, """
  <script type="text/ecmascript">
    <![CDATA[

	function init(evt)
	{
	    if ( window.svgDocument == null )
	    {
		svgDocument = evt.target.ownerDocument;
	    }

	    tooltip = svgDocument.getElementById('tooltip');
	    tooltip_bg = svgDocument.getElementById('tooltip_bg');

	}

	function ShowTooltip(evt, mouseovertext)
	{
	    var position = document.rootElement.createSVGPoint();
	    position.x = evt.clientX;
	    position.y = evt.clientY;
	    var matrix = document.rootElement.getScreenCTM();
	    var correctPosition=position.matrixTransform(matrix.inverse());
	    var mousex = correctPosition.x;
	    var mousey = correctPosition.y;
	    tooltip.setAttributeNS(null,"x",mousex+4);
	    tooltip.setAttributeNS(null,"y",mousey+14);
	    tooltip.firstChild.data = mouseovertext;
	    tooltip.setAttributeNS(null,"visibility","visible");

	    length = tooltip.getComputedTextLength();
	    tooltip_bg.setAttributeNS(null,"width",length+8);
	    tooltip_bg.setAttributeNS(null,"x",mousex+1);
	    tooltip_bg.setAttributeNS(null,"y",mousey+1);
	    tooltip_bg.setAttributeNS(null,"visibility","visible");
	}

	function HideTooltip(evt)
	{
	    tooltip.setAttributeNS(null,"visibility","hidden");
	    tooltip_bg.setAttributeNS(null,"visibility","hidden");
	}

    ]]>
  </script>
  """)

  for i=1:numpoints
    j = paired[i]
    if j > i && (length(visible) == 0 || (i in visible && j in visible))
      #write(svgout, string("<line id=\"bond_", i,"\" x1=\"", points[i][1], "\" y1=\"",points[i][2],"\"  x2=\"", points[j][1],"\" y2=\"", points[j][2] ,"\" style=\"stroke-width:", 8.0,";stroke:", rgbstring(getcolor(linkgradient, linkvalues[i])),"\""))

      dist = sqrt((points[i][1]-points[j][1])^2.0 + (points[i][2]-points[j][2])^2.0)
      linewidth = 16.0
      if dist > 60.0
        linewidth = 8.0
      end
      #println("DIST=",dist)
      write(svgout, string("<path id=\"bond_", i,"\" d=\"M", points[i][1], ",",points[i][2]," L", points[j][1],",", points[j][2] ,"\" style=\"fill:none;stroke-width:", linewidth,";stroke:", rgbstring(getcolor(linkgradient, linkvalues[i]/maxcoevolutionvalue)),"\""))
      #write(svgout, "/>\n")
      #write(svgout, string(""" cursor="help" onmousemove="ShowTooltip(evt, '""", @sprintf("Posterior probability γ≠0: %0.3f", linkvalues[i]),";  ",@sprintf("posterior probability paired: %0.3f",linkvalues2[i]),"""')" onmouseout="HideTooltip(evt)"/>\n"""))
      write(svgout, string(""" cursor="help" onmousemove="ShowTooltip(evt, '""", @sprintf("Coevolution: %0.3f", linkvalues[i]),"""')" onmouseout="HideTooltip(evt)"/>\n"""))
    end
  end

  for i=1:numpoints-1
    if (length(visible) == 0 || (i in visible && (i+1) in visible))
      x1 = points[i][1]
      y1 = points[i][2]
      x2 = points[i+1][1]
      y2 = points[i+1][2]
      write(svgout, string("<path id=\"covalent", i, "_", i+1,"\" d=\"M", x1, " ", y1, " L ", x2," ", y2 ,"\" style=\"stroke-width:", 6.0,";stroke:lightgray;fill:none;marker-mid:url(#markerArrow)\"/>\n"))
    end
  end

  minbackvalue = minimum(backvalues)
  maxbackvalue = maximum(backvalues)
  for i=1:numpoints
    if (length(visible) == 0 || (i in visible))
      x = points[i][1]
      y = points[i][2]
      nucbackgroundcolor = "#dd3333"
      nucbackgroundcolor = "#ffffff"
      if !isnan(backvalues[i])
        denominator = maxbackvalue-minbackvalue
        #v = (backvalues[i]-minbackvalue)/(denominator == 0.0 ? 1.0 : denominator)
        v = min(max(0.0, backvalues[i]),1.0)
        nucbackgroundcolor = rgbstring(getcolor(backgradient, v))
      end

      if i == 1 || i == numpoints || (startpos+i-1) % 5 == 0
          if i <= div(numpoints, 2)
              write(svgout, "<path d=\"M$(x-radius-3.0) $(y) L$(x-radius+3.0) $(y)\" style=\"stroke:black;stroke-width:5.0;\"/>")
              #write(svgout, "<text x=\"$(x-radius-3.5)\" y=\"$(y+2.0)\" alignment-baseline=\"middle\" style=\"font-size:18px\" text-anchor=\"end\">$(startpos+i-1)</text>")
              write(svgout, "<text x=\"$(x-radius-3.5)\" y=\"$(y+radius-9.0)\" style=\"font-size:18px\" text-anchor=\"end\">$(startpos+i-1)</text>")

          else
              write(svgout, "<path d=\"M$(x+radius-3.0) $(y) L$(x+radius+3.0) $(y)\" style=\"stroke:black;stroke-width:5.0;\"/>")
              #write(svgout, "<text x=\"$(x+radius+3.5)\" y=\"$(y+2.0)\" alignment-baseline=\"middle\" style=\"font-size:18px\" text-anchor=\"start\">$(startpos+i-1)</text>")
              write(svgout, "<text x=\"$(x+radius+3.5)\" y=\"$(y+radius-9.0)\" style=\"font-size:18px\" text-anchor=\"start\">$(startpos+i-1)</text>")
          end
      end

      if gapfrequency[i] > 0.5
        #<rect x="50" y="20" rx="20" ry="20" width="150" height="150"
        write(svgout, string("<rect id=\"nucleotide_",i,"\" x=\"",x-radius,"\" y=\"",y-radius,"\" width=\"",radius*2.0,"\" height=\"",radius*2.0,"\" rx=\"10.0\" ry=\"10.0\" stroke-dasharray=\"5, 5\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
        #write(svgout, string("<circle id=\"nucleotide_",i,"\" cx=\"",x,"\" cy=\"",y,"\" r=\"",radius,"\" stroke-dasharray=\"5, 5\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
      else
        write(svgout, string("<rect id=\"nucleotide_",i,"\" x=\"",x-radius,"\" y=\"",y-radius,"\" width=\"",radius*2.0,"\" height=\"",radius*2.0,"\" rx=\"10.0\" ry=\"10.0\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
        #write(svgout, string("<circle id=\"nucleotide_",i,"\" cx=\"",x,"\" cy=\"",y,"\" r=\"",radius,"\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
      end
      #write(svgout, "/>\n")
      write(svgout, string(""" cursor="help" onmousemove="ShowTooltip(evt, '""","Position ",startpos+i-1," (",get(mapping,i+startpos-1,"NA"),");  $(backlabel): ",@sprintf("%0.2f; gaps: %0.1f%%",backvalues[i], gapfrequency[i]*100.0),"""')" onmouseout="HideTooltip(evt)"/>\n"""))
    end
  end


  for i=1:size(data,2)
    if (length(visible) == 0 || (i in visible))
      v = zeros(Float64,4)
      for j=1:size(data,1)
        w = data[j,i,:]/sum(data[j,i,:])
        for k=1:4
          v[k] += w[k]
        end
      end
      v /= sum(v)
      h = zeros(Float64,4)
      for k=1:4
        h[k] -= v[k]*log(v[k])/log(2.0)
      end
      e = (3.0/log(2.0))/size(data,2)
      for k=1:4
        v[k] = v[k]*(2.0 - h[k] - e)
        if isnan(v[k])
          v[k] = 0.0
        end
      end
      v /= 2.0

      #println(h,"\t",v)
      x = points[i][1]
      y = points[i][2]
      logo_offset = 1.0*scale
      logoyoffset = (1.0 - sum(v))*radius
      additionaltags = string(""" cursor="help" onmousemove="ShowTooltip(evt, '""","Position ",startpos+i-1," (",get(mapping,i+startpos-1,"NA"),"); $(backlabel): ",@sprintf("%0.2f; gaps: %0.1f%%",backvalues[i], gapfrequency[i]*100.0),"')\" onmouseout=\"HideTooltip(evt)\"")
      #additionaltags = ""

      uselight = false
      if !isnan(backvalues[i])
        backcolor = getcolor(backgradient, min(max(0.0, backvalues[i]),1.0))
        backaverage = (backcolor[1] + backcolor[2] + backcolor[3])/3.0
        if backaverage <= 127.0
          uselight = true
        end
        if isnan(backvalues[i])
          uselight = false
        end
      end
      write(svgout, string(sequencelogo(string("logo",i), x, y-radius+logo_offset+logoyoffset, radius*2.0, radius*2.0-logo_offset*3.0, v, additionaltags, uselight),"\n"))
    end
  end

  write(svgout,"""
  <rect class="tooltip_bg" id="tooltip_bg"
      x="0" y="0" rx="4" ry="4"
      width="55" height="17" visibility="hidden"/>
  <text class="tooltip" id="tooltip"
      x="0" y="0" visibility="hidden">Tooltip</text>
  """)

  write(svgout, "</svg>\n")
  return points
end



function loadReplicationCapacityData(dataset, alignmentfile, infile)
  #sequence, pairedsites = readctfile("/media/michael/Sandisk500GB/data/hiv1-SHAPE-MaP.ct")
  #sequence = "GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGAAAGTAAAGCCAGAGGAGATCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCGGTATTAAGCGGGGGAGAATTAGATAAATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAACAATATAAACTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTTTTAGAGACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAATAGCAGTCCTCTATTGTGTGCATCAAAGGATAGATGTAAAAGACACCAAGGAAGCCTTAGATAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAGGCACAGCAAGCAGCAGCTGACACAGGAAACAACAGCCAGGTCAGCCAAAATTACCCTATAGTGCAGAACCTCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTAATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAATACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGATTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACACATAATCCACCTATCCCAGTAGGAGAAATCTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAAGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCATGTCAGGGAGTGGGGGGACCCGGCCATAAAGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATCCAGCTACCATAATGATACAGAAAGGCAATTTTAGGAACCAAAGAAAGACTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACATAGCCAAAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCCACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAAGAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAGCTTCCCTCAGATCACTCTTTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGCGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGCTGCACTTTAAATTTTCCCATTAGTCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGAAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGATTTCTGGGAAGTTCAATTAGGAATACCACATCCTGCAGGGTTAAAACAGAAAAAATCAGTAACAGTACTGGATGTGGGCGATGCATATTTTTCAGTTCCCTTAGATAAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAGTGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTCATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAGGACAGCTGGACTGTCAATGACATACAGAAATTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTATGCAGGGATTAAAGTAAGGCAATTATGTAAACTTCTTAGGGGAACCAAAGCACTAACAGAAGTAGTACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGGGAGATTCTAAAAGAACCGGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAGTATGCAAGAATGAAGGGTGCCCACACTAATGATGTGAAACAATTAACAGAGGCAGTACAAAAAATAGCCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAATTACCCATACAAAAGGAAACATGGGAAGCATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTCAATACCCCTCCCTTAGTGAAGTTATGGTACCAGTTAGAGAAAGAACCCATAATAGGAGCAGAAACTTTCTATGTAGATGGGGCAGCCAATAGGGAAACTAAATTAGGAAAAGCAGGATATGTAACTGACAGAGGAAGACAAAAAGTTGTCCCCCTAACGGACACAACAAATCAGAAGACTGAGTTACAAGCAATTCATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTGACAGACTCACAATATGCATTGGGAATCATTCAAGCACAACCAGATAAGAGTGAATCAGAGTTAGTCAGTCAAATAATAGAGCAGTTAATAAAAAAGGAAAAAGTCTACCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTGGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGAAGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTACCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGGGAAGCCATGCATGGACAAGTAGACTGTAGCCCAGGAATATGGCAGCTAGATTGTACACATTTAGAAGGAAAAGTTATCTTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTAATTCCAGCAGAGACAGGGCAAGAAACAGCATACTTCCTCTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAGTACATACAGACAATGGCAGCAATTTCACCAGTACTACAGTTAAGGCCGCCTGTTGGTGGGCGGGGATCAAGCAGGAATTTGGCATTCCCTACAATCCCCAAAGTCAAGGAGTAATAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAGATCCAGTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATCAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAACACATGGAAAAGATTAGTAAAACACCATATGTATATTTCAAGGAAAGCTAAGGACTGGTTTTATAGACATCACTATGAAAGTACTAATCCAAAAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAAATTAGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGACCTAGCAGACCAACTAATTCATCTGCACTATTTTGATTGTTTTTCAGAATCTGCTATAAGAAATACCATATTAGGACGTATAGTTAGTCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAGTACTTGGCACTAGCAGCATTAATAAAACCAAAACAGATAAAGCCACCTTTGCCTAGTGTTAGGAAACTGACAGAGGACAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCATACAATGAATGGACACTAGAGCTTTTAGAGGAACTTAAGAGTGAAGCTGTTAGACATTTTCCTAGGATATGGCTCCATAACTTAGGACAACATATCTATGAAACTTACGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATGACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAGTAAGTAGTACATGTAATGCAACCTATAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAGTATCAGCACTTGTGGAGATGGGGGTGGAAATGGGGCACCATGCTCCTTGGGATATTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGATAAGGTGCAGAAAGAATATGCATTCTTTTATAAACTTGATATAGTACCAATAGATAATACCAGCTATAGGTTGATAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATCAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGATGTAGTAATTAGATCTGCCAATTTCACAGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGTATCCGTATCCAGAGGGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATGCCACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACACTCCCATGCAGAATAAAACAATTTATAAACATGTGGCAGGAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACTGGGCTGCTATTAACAAGAGATGGTGGTAATAACAACAATGGGTCCGAGATCTTCAGACCTGGAGGAGGCGATATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGATATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAACAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATAACATGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAATCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTAGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAACTTGCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTATTACAAGCAGCTTATAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTGCTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATGGGGTGGGAGCAGTATCTCGAGACCTAGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTAACAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAAGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGGGGTCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGGTAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCTGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATGCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAA"
  #=
  mapping, revmapping = createmapping2("/media/michael/Sandisk500GB/data/NL4-3.fasta", "/media/michael/Sandisk500GB/data/pNL4-3_trunc.fasta")
  for i=1:length(mapping)
    println(i,"\t",mapping[i])
  end=#
  pairs = ["GC", "CG", "AT", "TA", "GT", "TG"]
  sequence, pairedsites = readctfile("/media/michael/Sandisk500GB/data/hiv1-SHAPE-MaP.ct")
  mapping, revmapping = createmapping(alignmentfile, sequence)

  table = readtable(infile)
  len = length(table[1])
  data = Dict{Int,Array{Float64,1}}()
  for i=1:len
    name = table[1][i]
    nuc1 = table[2][i][1]
    nuc2 = table[2][i][end]
    pos = parse(Int,table[2][i][2:end-1])
    rc = table[6][i]
    j = pos - 454
    nuc3 = 0
    nuc4 = 0
    nuc4 = sequence
    if j > 0 && pairedsites[j] != 0
      pair = string(sequence[j],sequence[pairedsites[j]])
      mutpair = string(nuc2,sequence[pairedsites[j]])
      breakspair = !(mutpair in pairs)
      #breakspair = true
      # && !(name != "5’LTR" && name != "3’LTR")
      if breakspair
        index = get(revmapping, j, 0)
        if pairedsites[j] < j && pairedsites[j] > 0
          index = get(revmapping, pairedsites[j], 0)
        end
        println(name, "\t", pos,"\t",j,"\t",nuc1,"->",nuc2, "\t", pair,"\t",mutpair,"\t",breakspair,"\t", rc)
        #index = min(i, pairedsites[i])
        if index > 0
          cur = get(data,index,[])
          push!(cur,rc)
          data[index] = cur
        end
      end
    end
    #println(nuc1,"\t",nuc2,"\t",pos)
  end

  res = zeros(Float64,dataset.numcols)
  for i=1:dataset.numcols
    if i in keys(data)
      res[i] = mean(data[i])
    else
      res[i] = NaN
    end
  end

  #println(res)
  return res
end

function loadUnmappedNL43ShapeReactivities2()
  sequence, paired = readctfile("/media/michael/Sandisk500GB/data/hiv1-SHAPE-MaP.ct")
  shapereactivities = zeros(Float64, length(sequence))
  shapeprobs = zeros(Float64, length(sequence))
  table = readtable("/media/michael/Sandisk500GB/data/NL4-3 SHAPE reactivities.csv")

  vec1 = table[3]
  vec2 = table[4]
  for j=1:length(sequence)
    k = paired[j]
    if j > 0  && k > 0 && j <= length(vec1) && k <= length(vec1) && !isna(vec1[j]) && !isna(vec1[k])
        shapereactivities[j]  =  (vec1[j]+vec1[k])/2.0
    else
      shapereactivities[j] = NaN
    end

    if j > 0  && k > 0 && j <= length(vec2) && k <= length(vec2) && !isna(vec2[j]) && !isna(vec2[k])
      shapeprobs[j]  =  (vec2[j]*vec2[k])
    else
      shapeprobs[j] = NaN
    end
  end

  return shapereactivities,shapeprobs
end

function loadUnmappedNL43ShapeReactivities()
  sequence, paired = readctfile("/media/michael/Sandisk500GB/data/hiv1-SHAPE-MaP.ct")
  shapereactivities = zeros(Float64, length(sequence))
  shapeprobs = zeros(Float64, length(sequence))
  table = readtable("/media/michael/Sandisk500GB/data/NL4-3 SHAPE reactivities.csv")

  vec1 = table[3]
  vec2 = table[4]
  for j=1:length(sequence)
    if isna(vec1[j])
      shapereactivities[j] = NaN
    else
      shapereactivities[j] = vec1[j]
    end

    if isna(vec2[j])
      shapeprobs[j] = NaN
    else
      shapeprobs[j] = vec2[j]
    end
  end

  return shapereactivities,shapeprobs
end

function loadNL43ShapeReactivities(dataset, alignmentfile, sequence)
  table = readtable("/media/michael/Sandisk500GB/data/NL4-3 SHAPE reactivities.csv")
  shapereactivities = zeros(Float64, dataset.numcols)
  shapeprobs = zeros(Float64, dataset.numcols)
  vec1 = table[3]
  #vec1 = loadReplicationCapacityData("/media/michael/Sandisk500GB/data/12977_2014_124_MOESM6_ESM.csv")
  vec2 = table[4]
  sequence2, pairedsites2 = readctfile("/media/michael/Sandisk500GB/data/hiv1-SHAPE-MaP.ct")
  mapping, revmapping = createmapping(alignmentfile, sequence2)
  for i=1:dataset.numcols
    j = get(mapping, i, 0)
    k = 0
    if j > 0 && pairedsites2[j] > j
      k = pairedsites2[j]
    end
    if j > 0  && k > 0 && j <= length(vec1) && k <= length(vec1) && !isna(vec1[j]) && !isna(vec1[k])
        shapereactivities[i]  =  (vec1[j]+vec1[k])/2.0
    else
      shapereactivities[i] = NaN
    end

    if j > 0  && k > 0 && j <= length(vec2) && k <= length(vec2) && !isna(vec2[j]) && !isna(vec2[k])
      shapeprobs[i]  =  (vec2[j]*vec2[k])
    else
      shapeprobs[i] = NaN
    end
  end
  return shapereactivities,shapeprobs,mapping, revmapping
end

function loadNL43ShapeReactivities2(dataset, alignmentfile, ctFile, csvfile)
  sequence, dummy = readctfile(ctFile)
  table = CSV.read(csvfile)
  shapereactivities = zeros(Float64, dataset.numcols)
  shapeprobs = zeros(Float64, dataset.numcols)
  vec1 = table[3]
  #vec1 = loadReplicationCapacityData("/media/michael/Sandisk500GB/data/12977_2014_124_MOESM6_ESM.csv")
  vec2 = table[4]
  sequence2, pairedsites2 = readctfile(ctFile)
  mapping, revmapping = createmapping(alignmentfile, sequence2)
  for i=1:dataset.numcols
    j = get(mapping, i, 0)
    if j > 0 && j <= length(vec1) && !ismissing(vec1[j])
        shapereactivities[i]  =  vec1[j]
    else
      shapereactivities[i] = NaN
    end

    if j > 0 && j <= length(vec2) && !ismissing(vec2[j])
      shapeprobs[i]  = vec2[j]
    else
      shapeprobs[i] = NaN
    end
  end
  return shapereactivities,shapeprobs,mapping, revmapping
end

function loadHCVShapeReactivities(dataset, alignmentfile, sequence)
  table = readtable("/media/michael/Sandisk500GB/data/hcvshape/con1.csv")
  shapereactivities = zeros(Float64, dataset.numcols)
  shapeprobs = zeros(Float64, dataset.numcols)
  vec1 = table[3]
  sequence2, pairedsites2 = readctfile("/media/michael/Sandisk500GB/data/hcvshape/HCV_Con1b_genome_model.ct")
  mapping, revmapping = createmapping(alignmentfile, sequence2)
  for i=1:dataset.numcols
    j = get(mapping, i, 0)
    k = 0
    if j > 0 && pairedsites2[j] > j
      k = pairedsites2[j]
    end
    #if j > 0  && j <= length(vec1)&& !isna(vec1[j]) && vec1[j] != -999.0
      #shapereactivities[i]  =  vec1[j]
    if j > 0  && k > 0 && j <= length(vec1) && k <= length(vec1) && !isna(vec1[j]) && !isna(vec1[k]) && vec1[j] != -999.0  && vec1[k] != -999.0
      shapereactivities[i]  =  (vec1[j]+vec1[k])/2.0
      println(i,"\t",j,"\t",sequence2[j],sequence2[k],"\t",shapereactivities[i])
    else
      shapereactivities[i] = NaN
    end
  end
  return shapereactivities,shapereactivities
end

function loadShapeReactivities(dataset, alignmentfile, sequence2, pairedsites2, csvfile, col::Int=3)
  table = readtable(csvfile)
  shapereactivities = zeros(Float64, dataset.numcols)
  allreactivities = zeros(Float64, dataset.numcols)
  vec1 = table[col]
  mapping, revmapping = createmapping(alignmentfile, sequence2)

  for i=1:dataset.numcols
    j = get(mapping, i, 0)
    k = 0
    if j > 0 && pairedsites2[j] > j
      k = pairedsites2[j]
    end

    if j > 0 && !isna(vec1[j])
      allreactivities[i] = vec1[j]
    else
      allreactivities[i] = NaN
    end

    if j > 0  && k > 0 && j <= length(vec1) && k <= length(vec1) && !isna(vec1[j]) && !isna(vec1[k]) && vec1[j] != -999.0  && vec1[k] != -999.0
      shapereactivities[i]  =  (vec1[j]+vec1[k])/2.0
      println(i,"\t",j,"\t",sequence2[j],sequence2[k],"\t",shapereactivities[i])
    else
      shapereactivities[i] = NaN
    end
  end
  return shapereactivities, allreactivities
end
