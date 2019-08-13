#include("SubstructureRanking.jl")
using Printf
include("Viridis.jl")

function pairedsitesgraph(pairedsites::Array{Int,1})
  fout = open("outfile.json", "w")
  write(fout, "{\n")
  write(fout, "\"nodes\": [")
  for i=1:length(pairedsites)
    write(fout, string("{\"id\": \"pos", i,"\", \"group\": \"1\"},\n"))
  end
  write(fout, "],\n")
  write(fout, "\"links\": [")
  for i=1:length(pairedsites)
    if i < length(pairedsites)
      write(fout, string("{\"source\": \"pos", i,"\", \"target\": \"pos",i+1,"\", \"value\": \"2\"},\n"))
    end
    if pairedsites[i] > i
      write(fout, string("{\"source\": \"pos", i,"\", \"target\": \"pos",pairedsites[i],"\", \"value\": \"1\"},\n"))
    end
  end
  write(fout, "]\n\n")
  write(fout, "}\n")
  close(fout)
end

function sequencelogo(id::AbstractString, x::Float64, y::Float64, width::Float64, height::Float64, h2::Array{Float64,1}, additionaltags::AbstractString="", uselight::Bool=false)
  h = h2 / sum(h2[1:4])

  fontHeight = 64.0
  scale = height / fontHeight
  base = y
  rna = "ACGU"
  ret = ""

  imagesdark = ["A3.svg", "C3.svg", "G3.svg", "U3.svg"]
  images = imagesdark
  if uselight
    images = ["A3_light.svg", "C3_light.svg", "G3_light.svg", "U3_light.svg"]
  end

  for i=1:length(h)
    fontHeightScale = h[i]
    b = rna[i]
    #base += fontHeightScale * scale * fontHeight
    xf = ((x-12.0) / scale)
    yf = base / (fontHeightScale*scale)
    if xf < 1e9 && yf < 1e9
      #<image transform="scale(1.0,1.0)" xlink:href="A1.png" x="0" y="0"/>transform="
      ret = string(ret, "<image ",additionaltags," transform=\"scale(", scale, ",",fontHeightScale * scale * 1.0, ")\" x=\"", xf, "\" y=\"", (yf / 1.0), "\" id=\"", id, "_", i, "\"   xlink:href=\"",images[i],"\">")
      #ret = string(ret, "<tspan id=\"""base""\">", b, "</tspan>\n")
    #  ret = string(ret, "<tspan>", b, "</tspan>")
      ret = string(ret, "</image>\n")
    end
    base += fontHeightScale * scale * fontHeight
  end

  return ret
end

#=
function sequencelogo(id::AbstractString, x::Float64, y::Float64, width::Float64, height::Float64, h::Array{Float64,1}, additionaltags::AbstractString="")
  fontHeight = 32.0
  scale = height / fontHeight
  base = y
  rna = "ACGU"
  ret = ""
  for i=1:length(h)
    fontHeightScale = h[i]
    b = rna[i]
    base += fontHeightScale * scale * fontHeight
    xf = x / scale
    yf = base / (fontHeightScale*scale)
    if xf < 128000 && yf < 128000
      ret = string(ret, "<text ",additionaltags," transform=\"scale(", scale, ",",fontHeightScale * scale * 1.1, ")\" x=\"", xf, "\" y=\"", (yf / 1.1), "\" id=\"", id, "_", i, "\"  style=\"font-size:", floor(Int, fontHeight), ";stroke:none;fill:", "black",";text-anchor:middle;font-family: Helvetica, Arial, Sans-Serif;\">")
      #ret = string(ret, "<tspan id=\"""base""\">", b, "</tspan>\n")
      ret = string(ret, "<tspan>", b, "</tspan>")
      ret = string(ret, "</text>\n")
    end

  end

  return ret
end=#

#=
function generatesubstructures(paired::Array{Int,1}, minlength::Int, maxsublength::Int, maxlength::Int)
  endpoints = Tuple{Int,Int}[]
  len = length(paired)
  i = 1
  startindex = 1
  endindex = 1
  while i <= len
    if paired[i] > i
      sublen = paired[i]-i+1
      if sublen <= maxsublength
        startindex = i
        endindex = paired[i]
        push!(endpoints, (startindex,endindex))
        i = endindex
      end
    end
    i += 1
  end
  endindex = endpoints[end][end]+1
  println("Z",paired[endindex:end], "\t", endpoints)

  finalendpoints = Tuple{Int,Int}[]
  startr = 1
  endr = 1
  r = 1
  while r <= length(endpoints)
    sublen = endpoints[r][2] - endpoints[startr][1] + 1
    println(endpoints)
    println(r,"\t", startr, "\t", endr, "\t", sublen, "\t", maxlength)
    if sublen <= maxlength

    else
      push!(finalendpoints, (endpoints[startr][1], endpoints[r-1][2]))
      startr = r
    end
    r += 1
  end
  if startr <= length(endpoints)
    push!(finalendpoints, (endpoints[startr][1], endpoints[end][2]))
  end
  println("Y", finalendpoints, "\t", startr, "\t", length(endpoints))

  substructures = Array{Int,1}[]
  starta = 1
  enda = 1
  finalendpoints2 = Tuple{Int,Int}[]
  for (startstruct endindex) in finalendpoints
    enda = startindex-1
    zlen = enda-starta+1
    if zlen > 0
      push!(finalendpoints2, (starta, enda))
      push!(substructures, zeros(Int,zlen))
      println(starta,"\t",enda,"\t", getdotbracketstring(substructures[end]))
    end
    substructure = copy(paired[startindex:endindex])
    for j=1:length(substructure)
      if substructure[j] != 0
        substructure[j] = substructure[j] - startindex + 1
      elseif !(startindex <= substructure[j] <= endindex)
        substructure[j]  = 0
      end
    end
    starta = endindex+1

    push!(finalendpoints2, (startindex, endindex))
    println(startindex,"\t",endindex,"\t", getdotbracketstring(substructure))
    push!(substructures, substructure)
  end
  starta = finalendpoints[end][2]+1
  enda = length(paired)
  zlen = enda-starta+1
  if enda-starta > 0
    push!(finalendpoints2, (starta, enda))
    push!(substructures, zeros(Int,zlen))
    println(starta,"\t",enda,"\t", getdotbracketstring(substructures[end]))
  end
  println("A", length(finalendpoints2), "\t", length(substructures))

  for substructure in substructures
    for i=1:length(substructure)
      substructure[i] = max(0, substructure[i])
    end
  end

  return finalendpoints2, substructures
end=#

function get2Dpoints(paired::Array{Int,1})

  temp = Array{Float64,1}[]
  allzero = true
  for i=1:length(paired)
    push!(temp, Float64[i*18.0,0.0])
    if paired[i] != 0
      allzero = false
      break
    end
  end
  if allzero
    return temp, 18.0, 0.0, 18.0*length(paired),0.0
  end
  dbnstring = getdotbracketstring(paired)
  binarypath = abspath(joinpath(@__DIR__, "../binaries/RNAVisualisation.jar"))
  data = read(`java -jar $(binarypath) "$dbnstring"`, String)
  ret = split(data,"\n")
  points = Array{Float64,1}[]
  minx = Inf
  miny = Inf
  maxx = -Inf
  maxy = -Inf
  for line in ret
    spl = split(strip(line),",")
    if length(spl) == 2
      p1 = parse(Float64,spl[1])
      p2 = parse(Float64,spl[2])
      minx = min(minx,p1)
      miny = min(miny,p2)
      maxx = max(maxx,p1)
      maxy = max(maxy,p2)
      push!(points, Float64[p1, p2])
    end
  end
  dy = points[1][2]
  for i=1:length(points)
    points[i] -= Float64[minx,dy]
  end
  return points, 0.0, miny-dy, maxx-minx, maxy-dy
end

function allunpaired(paired::Array{Int,1})
  for i=1:length(paired)
    if paired[i] != 0
      return false
    end
  end
  return true
end

function drawunpaired(len::Int, k::Int)
  scale = 1.0
  radius = 7.5*scale

  x = 0.0
  y = 0.0

  points = Array{Float64}[]

  #push!(points, [startx, starty])
  pathstring = string("<path d=\"M ",x," ",y)
  i = 1
  m = 1.0
  printL = true
  for i=1:len
    push!(points, Float64[x,y])
    if printL
      pathstring = string(pathstring, " L")
    end
    pathstring = string(pathstring, " ",x," ",y)
    if i % k == 0
      pathstring = string(pathstring, " Q ",x+radius*1.5," ",y+m*radius*2.0)
      printL = false
      x += radius*3.0
      m *= -1.0
    else
      printL = true
      y += m*radius*2.5
    end
  end

  minx = Inf
  miny = Inf
  maxx = -Inf
  maxy = -Inf
  for p in points
    minx = min(minx,p[1])
    miny = min(miny,p[2])
    maxx = max(maxx,p[1])
    maxy = max(maxy,p[2])
  end
  return points,minx,miny,maxx,maxy
end

function splaystructure(paired::Array{Int,1}, minlength::Int, maxsublength::Int, maxlength::Int, visible::Array{Int,1}=Int[])
  finalendpoints, substructures = generatesubstructures(paired, minlength, maxsublength, maxlength)
  finalpoints = Array{Float64,1}[]
  addx = 0.0
  finalminx = Inf
  finalminy = Inf
  finalmaxx = -Inf
  finalmaxy = -Inf
  for (endpoints,substructure) in zip(finalendpoints,substructures)
    ##println("ZZZZ",endpoints)
    #
      #=
      substructure = zeros()
      for z=endpoints[1]:endpoints[2]

      end=#

      points, minx, miny, maxx, maxy = get2Dpoints(substructure)
      for point in points
        push!(finalpoints, point + Float64[addx, 0.0])
      end
    if length(visible) == 0 || (endpoints[1] in visible && endpoints[2] in visible)
      addx += maxx + 40.0
      finalminx = min(finalminx, minx)
      finalminy = min(finalminy, miny)
      finalmaxx = max(finalmaxx, addx)
      finalmaxy = max(finalmaxy, maxy)
    end
  end
  #=
  while length(finalpoints) != length(paired)
    push!(finalpoints, Float64[0.0, 0.0])
  end=#
  return finalpoints, finalminx, finalminy, finalmaxx, finalmaxy
end

function drawlegend(outfile, gradient, tickslabels)
  width = 700
  height = 700
  offsety = 20.0

  svgout = open(outfile,"w")
  write(svgout,"""
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" onload="init(evt)" width="$width" height="$height">
  """
  )

legendheight = 200
for i=1:legendheight
  v = (i-1.0)/(legendheight-1.0)
  nucbackgroundcolor = rgbstring(getcolor(gradient, 1.0 - v))
  write(svgout, string("<rect x=\"10\" y=\"",offsety+i-1,"\" width=\"20\" height=\"2\" style=\"stroke: none; fill:",nucbackgroundcolor,";\"/>\n"))
end

for t=1:length(tickslabels)
  p = (t-1.0)/(length(tickslabels)-1.0)
  pos = (1.0-p)*legendheight
  write(svgout, string("<path id=\"tickmark\" d=\"M", 26.0, " ", offsety+pos+1.0, " L ", 34.0," ", offsety+pos+1.0,"\" style=\"stroke-width:", 0.5,";stroke:black;fill:none\"/>\n"))
  write(svgout, string("<text x=\"37\" y=\"",offsety+pos+6.0,"\" style=\"font-size: 19px;\">",tickslabels[t],"</text>\n"))
end


 write(svgout, "</svg>\n")
end

function drawsplayedstructure(outfile, dataset::Dataset, paired::Array{Int,1}, backvalues::Array{Float64,1}, linkvalues::Array{Float64,1}, linkvalues2::Array{Float64,1}, mapping, visible::Array{Int,1}=Int[])
  maxcoevolutionvalue = 10.0
  ylgnbu = reverse([[1.0,1.0,0.85],[0.93,0.98,0.71],[0.78,0.91,0.71],[0.50,0.80,0.73],[0.27,0.71,0.76],[0.15,0.56,0.74],[0.15,0.37,0.64],[0.16,0.22,0.56],[0.06,0.14,0.34]])

  shapecolors = [(0,0,0,1.0),(0,0,0,1.0),  (33,177,79,1.0),(33,177,79,1.0),  (250,200,36,1.0),(250,200,36,1.0), (244,23,28,1.0),(244,23,28,1.0),(244,23,28,1.0)]
  shapepositions = [0.0, 0.3, 0.3, 0.4, 0.4, 0.6, 0.6, 1.0]
  shapegradient = ColorGradient(shapecolors, shapepositions)



  #linkgradient =  viridisgradient(ylgnbu)
  linkgradient =  viridisgradient(_viridis_data)
  #backgradient =  viridisgradient(_viridis_data)
  backgradient =  shapegradient
  #drawlegend(string(outfile, ".substitutionrate.svg"), backgradient, ["low", "moderate", "high"])
  #drawlegend(string(outfile, ".coevolutionprobability.svg"), linkgradient, ["0.00", "0.25", "0.50", "0.75", "1.00"])
  #drawlegend(string(outfile, ".shape.svg"), shapegradient, ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])
  points, minx, miny, maxx, maxy = splaystructure(paired, 1, 350, 150, visible)
  numpoints = length(points)

  scale = 2.0
  radius = 7.5*scale
  xoffset = radius
  yoffset = radius
  width = (maxx-minx)*scale + xoffset*2.0 + 0.0
  height = (maxy-miny)*scale + yoffset*2.0 + 20.0
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
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" onload="init(evt)" width="$width" height="$height">
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
      write(svgout, string(""" cursor="help" onmousemove="ShowTooltip(evt, '""", @sprintf("Posterior probability γ≠0: %0.3f", linkvalues[i]),";  ",@sprintf("posterior probability paired: %0.3f",linkvalues2[i]),"""')" onmouseout="HideTooltip(evt)"/>\n"""))
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
      if dataset.gapfrequency[i] > 0.5
        #<rect x="50" y="20" rx="20" ry="20" width="150" height="150"
        write(svgout, string("<rect id=\"nucleotide_",i,"\" x=\"",x-radius,"\" y=\"",y-radius,"\" width=\"",radius*2.0,"\" height=\"",radius*2.0,"\" rx=\"10.0\" ry=\"10.0\" stroke-dasharray=\"5, 5\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
        #write(svgout, string("<circle id=\"nucleotide_",i,"\" cx=\"",x,"\" cy=\"",y,"\" r=\"",radius,"\" stroke-dasharray=\"5, 5\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
      else
        write(svgout, string("<rect id=\"nucleotide_",i,"\" x=\"",x-radius,"\" y=\"",y-radius,"\" width=\"",radius*2.0,"\" height=\"",radius*2.0,"\" rx=\"10.0\" ry=\"10.0\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
        #write(svgout, string("<circle id=\"nucleotide_",i,"\" cx=\"",x,"\" cy=\"",y,"\" r=\"",radius,"\" style=\"stroke-width:2;stroke:black;fill:",nucbackgroundcolor,"\" "))
      end
      #write(svgout, "/>\n")
      write(svgout, string(""" cursor="help" onmousemove="ShowTooltip(evt, '""","Position ",i," (",get(mapping,i,0),");  ",@sprintf("mean posterior site rate: %0.3f;  gaps: %0.1f%%",backvalues[i], dataset.gapfrequency[i]*100.0),"""')" onmouseout="HideTooltip(evt)"/>\n"""))
    end
  end


  for i=1:size(dataset.data,2)
    if (length(visible) == 0 || (i in visible))
      v = zeros(Float64,4)
      for j=1:size(dataset.data,1)
        w = dataset.data[j,i,:]/sum(dataset.data[j,i,:])
        for k=1:4
          v[k] += w[k]
        end
      end
      v /= sum(v)
      h = zeros(Float64,4)
      for k=1:4
        h[k] -= v[k]*log(v[k])/log(2.0)
      end
      e = (3.0/log(2.0))/size(dataset.data,2)
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
      additionaltags = string(""" cursor="help" onmousemove="ShowTooltip(evt, '""","Position ",i," (",get(mapping,i,0),");  ",@sprintf("mean posterior site rate: %0.3f;  gaps: %0.1f%%",backvalues[i], dataset.gapfrequency[i]*100.0),"')\" onmouseout=\"HideTooltip(evt)\"")
      #additionaltags = ""

      uselight = false
      backcolor = getcolor(backgradient, min(max(0.0, backvalues[i]),1.0))
      backaverage = (backcolor[1] + backcolor[2] + backcolor[3])/3.0
      if backaverage <= 127.0
        uselight = true
      end
      if isnan(backvalues[i])
        uselight = false
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
