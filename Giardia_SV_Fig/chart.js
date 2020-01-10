var svg = d3.select('svg'),
margin = {top:20, right: 20, bottom:75, left: 60},
gSpace = +svg.attr('width')/3,
chrgap = 20,
width = +svg.attr('width') - margin.left - margin.right,
legendSpace = 100,
height = +svg.attr('height') - margin.top - margin.bottom - legendSpace;

var ref = [1491578, 1519671, 1958359, 2762469, 4471437],
sv_types = ["Duplications", "Deletions", "Insertions", "Inversions", "Inverted Duplications"],
awb = [],
ags = [],
abeav = [],
svwb = [],
svbeav = [],
svgs = [],
sorted_wb = [],
sorted_beav = [],
sorted_gs = [];

// Chart height must accomodate chr 3 and 4 (4720828 bp) plus gap between them
var c1size = height * ref[0]/(ref[2] + ref[3]) - chrgap/2,
c2size = height * ref[1]/(ref[2] + ref[3]) - chrgap/2,
c3size = height * ref[2]/(ref[2] + ref[3]) - chrgap/2,
c4size = height * ref[3]/(ref[2] + ref[3]) - chrgap/2,
c5size = height * ref[4]/(ref[2] + ref[3]) - chrgap/2;

var chr1 = d3.scaleLinear().domain([0, ref[0]]).rangeRound([c1size, 0]),
chr2 = d3.scaleLinear().domain([0, ref[1]]).rangeRound([c2size, 0]),
chr3 = d3.scaleLinear().domain([0, ref[2]]).rangeRound([c3size, 0]),
chr4 = d3.scaleLinear().domain([0, ref[3]]).rangeRound([c4size, 0]),
chr5 = d3.scaleLinear().domain([0, ref[4]]).rangeRound([c5size, 0]),
x1 = d3.scaleBand().domain([0, 1, 2]).rangeRound([0, width/3]).padding(0.1),
x2 = d3.scaleBand().domain([0, 1, 2]).rangeRound([0, width/3]).padding(0.1),
x3 = d3.scaleBand().domain([0, 1, 2]).rangeRound([0, width/3]).padding(0.1),
x4 = d3.scaleBand().domain([0, 1, 2]).rangeRound([0, width/3]).padding(0.1),
x5 = d3.scaleBand().domain([0, 1, 2]).rangeRound([0, width/3]).padding(0.1);

var min = 0,
max = 4471437;

var g1 = svg.append("g").attr("transform","translate(" + margin.left + "," + margin.top + ")");
var legend = svg.append("g").attr("transform", "translate(" + margin.left + "," + (height + margin.top + margin.bottom) + ")");

d3.queue()
	.defer(d3.tsv, "WB_hybrid_106_0150015723312338_1dsmartx0_humex_chrsonly.tab")
	.defer(d3.tsv, "GS_hybrid_gs3-20-2019_22372244_1dsmartx0_humex_chrsonly.tab")
	.defer(d3.tsv, "Beaver_hybrid_107218_2309_1dsmartx0_humex_chrsonly.tab")
	.defer(d3.tsv, "WB_allontreads_hybrid_106_0150015723312338smart1dx0_SVs_withOptions.vcf")
	.defer(d3.tsv, "Beaver_allontreads_hybrid_107218_2309_1dsmartx0_SVs_withOptions.vcf")
	.defer(d3.tsv, "GS_allontreads_hybrid_iseq_22372244_1dsmartx0_SVs_withOptions.vcf")
	.await(initialize);

function initialize(error, wb_align, gs_align, beaver_align, wb_sv, beav_sv, gs_sv) {
	if (error) {console.log(error); throw error;}
	
	setAxes();
	
	// Read data into appropriate arrays
	awb.push(wb_align);
	awb = awb[0];
	ags.push(gs_align);
	ags = ags[0];
	abeav.push(beaver_align);
	abeav = abeav[0];
	svwb.push(wb_sv);
	svwb = svwb[0];
	svbeav.push(beav_sv);
	svbeav = svbeav[0];
	svgs.push(gs_sv);
	svgs = svgs[0];
	
	//console.log(awb);
	//console.log(ags);
	//console.log(abeav);
	//console.log(svwb);
	
	// Sort each alignment data first by chromosome then by start position
	awb.sort(function(x,y){
		if (x.Chr > y.Chr){
			return 1;
		} else if (x.Chr < y.Chr) {
			return -1;
		}
		
		if (+x.Start < +y.Start) {
			return -1;
		} else if (+x.Start > +y.Start) {
			return 1;
		} else {
			return 0;
		}
	});
	ags.sort(function(x,y){
		if (x.Chr > y.Chr){
			return 1;
		} else if (x.Chr < y.Chr) {
			return -1;
		}
		
		if (+x.Start < +y.Start) {
			return -1;
		} else if (+x.Start > +y.Start) {
			return 1;
		} else {
			return 0;
		}
	});
	abeav.sort(function(x,y){
		if (x.Chr > y.Chr){
			return 1;
		} else if (x.Chr < y.Chr) {
			return -1;
		}
		
		if (+x.Start < +y.Start) {
			return -1;
		} else if (+x.Start > +y.Start) {
			return 1;
		} else {
			return 0;
		}
	});
	
	//console.log(awb);
	//console.log(ags);
	//console.log(abeav);
	
	mergeAdjacentRegions(awb, abeav, ags);
	//drawAlignments(sorted_wb, sorted_beav, sorted_gs);
	drawAlignments(awb, abeav, ags);
	drawLegend();
	//drawSVs(svwb, svbeav, svgs, sorted_wb, sorted_beav, sorted_gs);
	drawSVs(svwb, svbeav, svgs, awb, abeav, ags);
	styleVis();
	
}

/**
 * Create the axes used for mapping all alignment data
 */
function setAxes() {
	// Create each chromosome as its own y-axis positioned in 3 columns and create a corresponding x-axis
	g1.append("g").attr("transform", "translate(20," + height + ")").append("xaxis").call(d3.axisBottom(x1).tickSize(0));
	g1.append("g").attr("transform", "translate(20," + 0 + ")").call(d3.axisLeft(chr1).tickSize(0).tickValues([]));
	g1.append("g").attr("transform", "translate(0," + height + ")").append("xaxis").call(d3.axisBottom(x2).tickSize(0));
	g1.append("g").attr("transform", "translate(0," + (c1size + chrgap) + ")").call(d3.axisLeft(chr2).tickSize(0).tickValues([]));

	g1.append("g").attr("transform", "translate(" + gSpace + "," + height + ")").append("xaxis").call(d3.axisBottom(x3).tickSize(0));
	g1.append("g").attr("transform", "translate(" + gSpace + "," + 0 + ")").call(d3.axisLeft(chr3).tickSize(0).tickValues([]));
	g1.append("g").attr("transform", "translate(" + (gSpace - 20) + "," + height + ")").append("xaxis").call(d3.axisBottom(x4).tickSize(0));
	g1.append("g").attr("transform", "translate(" + (gSpace - 20) + "," + (c3size + chrgap) + ")").call(d3.axisLeft(chr4).tickSize(0).tickValues([]));
	
	g1.append("g").attr("transform", "translate(" + 2*gSpace + "," + height + ")").append("xaxis").call(d3.axisBottom(x5).tickSize(0));
	g1.append("g").attr("transform", "translate(" + 2*gSpace + "," + 0 + ")").call(d3.axisLeft(chr5).tickSize(0).tickValues([]));
	
	// Remove the line on each x-axis
	g1.selectAll("xaxis").remove();
	
	// Add labels to each chromosome
	g1.append("text").attr("x", -c1size/2).attr("y", -margin.left/4 + 20).attr("transform", "rotate(-90)")
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "middle").text("GLCHR01");
	g1.append("text").attr("x", -(c2size/2 + c1size + chrgap)).attr("y", -margin.left/4).attr("transform", "rotate(-90)")
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "middle").text("GLCHR02");
	g1.append("text").attr("x", -c3size/2).attr("y", -margin.left/4 + gSpace).attr("transform", "rotate(-90)")
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "middle").text("GLCHR03");
	g1.append("text").attr("x", -(c4size/2 + c3size + chrgap)).attr("y", -margin.left/4 + gSpace - 20).attr("transform", "rotate(-90)")
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "middle").text("GLCHR04");
	g1.append("text").attr("x", -c5size/2).attr("y", -margin.left/4 + 2*gSpace).attr("transform", "rotate(-90)")
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "middle").text("GLCHR05");
}

/**
 * Scans sorted alignment arrays to merge overlapping regions or regions separated by less than 1000 bp 
 * if they are part of the same query contig
 */
function mergeAdjacentRegions(wb, beav, gs) {
	sorted_wb.length = 0;
	sorted_beav.length = 0;
	sorted_gs.length = 0;
	
	var prevName = "";
	var firstRegion = "";
	
	for (var i = 0; i < wb.length; i++) {
		if (i == 0) {
			prevName = wb[i].QueryName;
			firstRegion = JSON.parse(JSON.stringify(wb[i]));
			continue;
		} else if (wb[i].QueryName == prevName) {
			// Same contig as previous, check if the regions overlap
			if (+wb[i].Start < (+wb[i-1].Start + +wb[i-1].Size)) {
				// Merge into previous region if it extends beyond current region
				if ((+wb[i].Start + +wb[i].Size) > (+firstRegion.Size + +firstRegion.Start)) {
					firstRegion.Size = String(+wb[i].Start + +wb[i].Size - (+firstRegion.Start));
				}
			} else if (+wb[i].Start < (+wb[i-1].Start + +wb[i-1].Size + 1001)) {
				// Space is < 1000, merge into previous region
				if ((+wb[i].Start + +wb[i].Size) > (+firstRegion.Size + +firstRegion.Start)) {
					firstRegion.Size = String(+wb[i].Start + +wb[i].Size - (+firstRegion.Start));
				}
			} else {
				// Space bigger, do not merge
				// Push previously built info
				if (+firstRegion.Size > 1000) {
					firstRegion.QuerySize = String(firstRegion.Size);
				}
				sorted_wb.push(firstRegion);
				
				// Load this new region in
				firstRegion = JSON.parse(JSON.stringify(wb[i]));
			}
		} else {
			if (+wb.Size <= 1000 || +wb.QuerySize <= 1000) {
				continue;
			}
			// Different contig, do not merge
			// Push previously built info
			if (+firstRegion.Size > 1000) {
				firstRegion.QuerySize = String(firstRegion.Size);
			}
			sorted_wb.push(firstRegion);

			// Reset first region to new contig
			firstRegion = JSON.parse(JSON.stringify(wb[i]));
		}
		
		prevName = wb[i].QueryName;
	}

	var prevName = "";
	var firstRegion = "";
	
	for (var i = 0; i < beav.length; i++) {
		if (i == 0) {
			prevName = beav[i].QueryName;
			firstRegion = JSON.parse(JSON.stringify(beav[i]));
			continue;
		} else if (beav[i].QueryName == prevName) {
			// Same contig as previous, check if the regions overlap
			if (+beav[i].Start < (+beav[i-1].Start + +beav[i-1].Size)) {
				// Merge into previous region if it extends beyond current region
				if ((+beav[i].Start + +beav[i].Size) > (+firstRegion.Size + +firstRegion.Start)) {
					firstRegion.Size = String(+beav[i].Start + +beav[i].Size - (+firstRegion.Start));
				}
			} else if (+beav[i].Start < (+beav[i-1].Start + +beav[i-1].Size + 1001)) {
				// Space is < 1000, merge into previous region
				if ((+beav[i].Start + +beav[i].Size) > (+firstRegion.Size + +firstRegion.Start)) {
					firstRegion.Size = String(+beav[i].Start + +beav[i].Size - (+firstRegion.Start));
				}
			} else {
				// Space bigger, do not merge
				// Push previously built info
				if (+firstRegion.Size > 1000) {
					firstRegion.QuerySize = String(firstRegion.Size);
				}
				sorted_beav.push(firstRegion);
				
				// Load this new region in
				firstRegion = JSON.parse(JSON.stringify(beav[i]));
			}
		} else {
			if (+beav.Size <= 1000 || +beav.QuerySize <= 1000) {
				continue;
			}
			// Different contig, do not merge
			// Push previously built info
			if (+firstRegion.Size > 1000) {
				firstRegion.QuerySize = String(firstRegion.Size);
			}
			sorted_beav.push(firstRegion);

			// Reset first region to new contig
			firstRegion = JSON.parse(JSON.stringify(beav[i]));
		}
		
		prevName = beav[i].QueryName;
	}
	
	var prevName = "";
	var firstRegion = "";
	
	for (var i = 0; i < gs.length; i++) {
		if (i == 0) {
			prevName = gs[i].QueryName;
			firstRegion = JSON.parse(JSON.stringify(gs[i]));
			continue;
		} else if (gs[i].QueryName == prevName) {
			// Same contig as previous, check if the regions overlap
			if (+gs[i].Start < (+gs[i-1].Start + +gs[i-1].Size)) {
				// Merge into previous region if it extends beyond current region
				if ((+gs[i].Start + +gs[i].Size) > (+firstRegion.Size + +firstRegion.Start)) {
					firstRegion.Size = String(+gs[i].Start + +gs[i].Size - (+firstRegion.Start));
				}
			} else if (+gs[i].Start < (+gs[i-1].Start + +gs[i-1].Size + 1001)) {
				// Space is < 1000, merge into previous region
				if ((+gs[i].Start + +gs[i].Size) > (+firstRegion.Size + +firstRegion.Start)) {
					firstRegion.Size = String(+gs[i].Start + +gs[i].Size - (+firstRegion.Start));
				}
			} else {
				// Space bigger, do not merge
				// Push previously built info
				if (+firstRegion.Size > 1000) {
					firstRegion.QuerySize = String(firstRegion.Size);
				}
				sorted_gs.push(firstRegion);
				
				// Load this new region in
				firstRegion = JSON.parse(JSON.stringify(gs[i]));
			}
		} else {
			if (+gs.Size <= 1000 || +gs.QuerySize <= 1000) {
				continue;
			}
			// Different contig, do not merge
			// Push previously built info
			if (+firstRegion.Size > 1000) {
				firstRegion.QuerySize = String(firstRegion.Size);
			}
			sorted_gs.push(firstRegion);

			// Reset first region to new contig
			firstRegion = JSON.parse(JSON.stringify(gs[i]));
		}
		
		prevName = gs[i].QueryName;
	}
	
	//console.log(sorted_wb);
	//console.log(sorted_beav);
	//console.log(sorted_gs);
}

 /**
  * Draw the genome alignments to the WB reference
  */
function drawAlignments(wb, beav, gs) {
	// AWB hybrid alignment
	g1.selectAll(".contig").data(wb.filter(function(d){return d.Chr == "GLCHR01" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x1(0) + 20;})
		.attr("y", function(d){return chr1(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr1(d.Start) - chr1(+d.Start + +d.Size);})
		.attr("class", "wb_contig")//.attr("num", function(d){return d.Num;})
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(wb.filter(function(d){return d.Chr == "GLCHR02" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x2(0);})
		.attr("y", function(d){return chr2(+d.Start + +d.Size) + c1size + chrgap;})
		.attr("width", 10)
		.attr("height", function(d){return chr2(d.Start) - chr2(+d.Start + +d.Size);})
		.attr("class", "wb_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(wb.filter(function(d){return d.Chr == "GLCHR03" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x3(0) + gSpace;})
		.attr("y", function(d){return chr3(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr3(d.Start) - chr3(+d.Start + +d.Size);})
		.attr("class", "wb_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(wb.filter(function(d){return d.Chr == "GLCHR04" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x4(0) + gSpace - 20;})
		.attr("y", function(d){return chr4(+d.Start + +d.Size) + c3size + chrgap;})
		.attr("width", 10)
		.attr("height", function(d){return chr4(d.Start) - chr4(+d.Start + +d.Size);})
		.attr("class", "wb_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(wb.filter(function(d){return d.Chr == "GLCHR05" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x5(0) + 2*gSpace;})
		.attr("y", function(d){return chr5(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr5(d.Start) - chr5(+d.Start + +d.Size);})
		.attr("class", "wb_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	
	// Beaver hybrid alignment
	g1.selectAll(".contig").data(beav.filter(function(d){return d.Chr == "GLCHR01" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x1(1) + 20;})
		.attr("y", function(d){return chr1(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr1(d.Start) - chr1(+d.Start + +d.Size);})
		.attr("class", "beaver_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(beav.filter(function(d){return d.Chr == "GLCHR02" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x2(1);})
		.attr("y", function(d){return chr2(+d.Start + +d.Size) + c1size + chrgap;})
		.attr("width", 10)
		.attr("height", function(d){return chr2(d.Start) - chr2(+d.Start + +d.Size);})
		.attr("class", "beaver_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(beav.filter(function(d){return d.Chr == "GLCHR03" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x3(1) + gSpace;})
		.attr("y", function(d){return chr3(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr3(d.Start) - chr3(+d.Start + +d.Size);})
		.attr("class", "beaver_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(beav.filter(function(d){return d.Chr == "GLCHR04" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x4(1) + gSpace - 20;})
		.attr("y", function(d){return chr4(+d.Start + +d.Size) + c3size + chrgap;})
		.attr("width", 10)
		.attr("height", function(d){return chr4(d.Start) - chr4(+d.Start + +d.Size);})
		.attr("class", "beaver_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(beav.filter(function(d){return d.Chr == "GLCHR05" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x5(1) + 2*gSpace;})
		.attr("y", function(d){return chr5(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr5(d.Start) - chr5(+d.Start + +d.Size);})
		.attr("class", "beaver_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	
	//BGS hybrid alignment
	g1.selectAll(".contig").data(gs.filter(function(d){return d.Chr == "GLCHR01" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x1(2) + 20;})
		.attr("y", function(d){return chr1(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr1(d.Start) - chr1(+d.Start + +d.Size);})
		.attr("class", "gs_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(gs.filter(function(d){return d.Chr == "GLCHR02" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x2(2);})
		.attr("y", function(d){return chr2(+d.Start + +d.Size) + c1size + chrgap;})
		.attr("width", 10)
		.attr("height", function(d){return chr2(d.Start) - chr2(+d.Start + +d.Size);})
		.attr("class", "gs_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(gs.filter(function(d){return d.Chr == "GLCHR03" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x3(2) + gSpace;})
		.attr("y", function(d){return chr3(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr3(d.Start) - chr3(+d.Start + +d.Size);})
		.attr("class", "gs_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(gs.filter(function(d){return d.Chr == "GLCHR04" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x4(2) + gSpace - 20;})
		.attr("y", function(d){return chr4(+d.Start + +d.Size) + c3size + chrgap;})
		.attr("width", 10)
		.attr("height", function(d){return chr4(d.Start) - chr4(+d.Start + +d.Size);})
		.attr("class", "gs_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
	g1.selectAll(".contig").data(gs.filter(function(d){return d.Chr == "GLCHR05" && d.Size > 1000 && d.QuerySize > 1000;})).enter().append("rect")
		.attr("x", function(d){return x5(2) + 2*gSpace;})
		.attr("y", function(d){return chr5(+d.Start + +d.Size);})
		.attr("width", 10)
		.attr("height", function(d){return chr5(d.Start) - chr5(+d.Start + +d.Size);})
		.attr("class", "gs_contig")
		.attr("start", function(d){return d.Start;}).attr("end", function(d){return +d.Start + +d.Size;})
		.attr("qstart", function(d){return d.QueryStart;}).attr("qend", function(d){return +d.QueryStart + +d.QuerySize;})
		.attr("id", function(d){return d.QueryName;}).attr("mapchr", function(d){return d.Chr;});
}

/**
 * Create the legend beneath the vis
 */
function drawLegend() {
	
	// Rectangles to coordinate colour
	legend.append("rect").attr("class", "wb_contig").attr("id", "legend")
		.attr("x", 0).attr("y", 0)
		.attr("width", 20).attr("height", 20);
	legend.append("rect").attr("class", "beaver_contig").attr("id", "legend")
		.attr("x", function(){return width/3;}).attr("y", 0)
		.attr("width", 20).attr("height", 20);
	legend.append("rect").attr("class", "gs_contig").attr("id", "legend")
		.attr("x", function(){return width*2/3;}).attr("y", 0)
		.attr("width", 20).attr("height", 20);
	
	// Labels
	legend.append("text").attr("x", 22).attr("y", 17)
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "start").text("AWB hybrid assembly");
	legend.append("text").attr("x", function(){return width/3 + 22;}).attr("y", 17)
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "start").text("beaver hybrid assembly");
	legend.append("text").attr("x", function(){return width*2/3 + 22;}).attr("y", 17)
		.style("fill","black").style("fill-opacity", "1")
		.style("font-size", "20").style("font-weight", "bold")
		.style("text-anchor", "start").text("BGS hybrid assembly");
}

/**
 * Add structural variants onto contigs
 */
function drawSVs(wb, beav, gs, dwb, dbeav, dgs) {
	// Add rectangles onto contig rectangles to represent the SVs
	// First find the contig rectangle that the SV should be drawn on
	for (var i = 0; i < wb.length; i++) {
		var mycontig = -1;
		var amt = 0;
		
		for (var j = 0; j < dwb.length; j++) {
			if ((dwb[j].QueryName !== wb[i].CHROM)) {
				continue;
			} else {				
				if (+dwb[j].QueryStart < +wb[i].POS) {					
					if ((+dwb[j].QuerySize + +dwb[j].QueryStart) > +wb[i].POS) {
						// This segment contains the start of the SV
						
						if ((+dwb[j].QuerySize + +dwb[j].QueryStart) > +wb[i].INFO.substring(wb[i].INFO.indexOf("END=") + 4, wb[i].INFO.indexOf(";STD_quant_start"))) {
							// This segment contains the end of the SV and so contains the whole SV
							mycontig = j;
							amt = 1;
						} else {
							// This SV is spread across multiple segments
							var myamt = ((+dwb[j].QuerySize + +dwb[j].QueryStart) - +wb[i].POS) / (+wb[i].INFO.substring(wb[i].INFO.indexOf("END=") + 4, wb[i].INFO.indexOf(";STD_quant_start")) - +wb[i].POS);
							if (myamt > amt) {
								mycontig = j;
								amt  = myamt;
							}
						}
					} else {
						// This segment is before the SV
						continue;
					}
				} else if ((+dwb[j].QuerySize + +dwb[j].QueryStart) < +wb[i].INFO.substring(wb[i].INFO.indexOf("END=") + 4, wb[i].INFO.indexOf(";STD_quant_start"))) {
					// This segment contains the end of the SV but not the beginning
					var myamt = (+wb[i].INFO.substring(wb[i].INFO.indexOf("END=") + 4, wb[i].INFO.indexOf(";STD_quant_start")) - +dwb[j].QueryStart) / (+wb[i].INFO.substring(wb[i].INFO.indexOf("END=") + 4, wb[i].INFO.indexOf(";STD_quant_start")) - +wb[i].POS);
					if (myamt > amt) {
						mycontig = j;
						amt  = myamt;
					}
				} else {
					// This segment is after the SV
					continue;
				}
			}
		} // end inner for
		
		// Store mycontig for drawing
		wb[i].contig = mycontig;		
	} // end outer for
	
	for (var i = 0; i < beav.length; i++) {
		var mycontig = -1;
		var amt = 0;
		
		for (var j = 0; j < dbeav.length; j++) {
			if ((dbeav[j].QueryName !== beav[i].CHROM)) {
				continue;
			} else {				
				if (+dbeav[j].QueryStart < +beav[i].POS) {					
					if ((+dbeav[j].QuerySize + +dbeav[j].QueryStart) > +beav[i].POS) {
						// This segment contains the start of the SV
						
						if ((+dbeav[j].QuerySize + +dbeav[j].QueryStart) > +beav[i].INFO.substring(beav[i].INFO.indexOf("END=") + 4, beav[i].INFO.indexOf(";STD_quant_start"))) {
							// This segment contains the end of the SV and so contains the whole SV
							mycontig = j;
							amt = 1;
						} else {
							// This SV is spread across multiple segments
							var myamt = ((+dbeav[j].QuerySize + +dbeav[j].QueryStart) - +beav[i].POS) / (+beav[i].INFO.substring(beav[i].INFO.indexOf("END=") + 4, beav[i].INFO.indexOf(";STD_quant_start")) - +beav[i].POS);
							if (myamt > amt) {
								mycontig = j;
								amt  = myamt;
							}
						}
					} else {
						// This segment is before the SV
						continue;
					}
				} else if ((+dbeav[j].QuerySize + +dbeav[j].QueryStart) < +beav[i].INFO.substring(beav[i].INFO.indexOf("END=") + 4, beav[i].INFO.indexOf(";STD_quant_start"))) {
					// This segment contains the end of the SV but not the beginning
					var myamt = (+beav[i].INFO.substring(beav[i].INFO.indexOf("END=") + 4, beav[i].INFO.indexOf(";STD_quant_start")) - +dbeav[j].QueryStart) / (+beav[i].INFO.substring(beav[i].INFO.indexOf("END=") + 4, beav[i].INFO.indexOf(";STD_quant_start")) - +beav[i].POS);
					if (myamt > amt) {
						mycontig = j;
						amt  = myamt;
					}
				} else {
					// This segment is after the SV
					continue;
				}
			}
		} // end inner for
		
		// Store mycontig for drawing
		beav[i].contig = mycontig;		
	} // end outer for
	
	for (var i = 0; i < gs.length; i++) {
		var mycontig = -1;
		var amt = 0;
		
		for (var j = 0; j < dgs.length; j++) {
			if ((dgs[j].QueryName !== gs[i].CHROM)) {
				continue;
			} else {				
				if (+dgs[j].QueryStart < +gs[i].POS) {					
					if ((+dgs[j].QuerySize + +dgs[j].QueryStart) > +gs[i].POS) {
						// This segment contains the start of the SV
						
						if ((+dgs[j].QuerySize + +dgs[j].QueryStart) > +gs[i].INFO.substring(gs[i].INFO.indexOf("END=") + 4, gs[i].INFO.indexOf(";STD_quant_start"))) {
							// This segment contains the end of the SV and so contains the whole SV
							mycontig = j;
							amt = 1;
						} else {
							// This SV is spread across multiple segments
							var myamt = ((+dgs[j].QuerySize + +dgs[j].QueryStart) - +gs[i].POS) / (+gs[i].INFO.substring(gs[i].INFO.indexOf("END=") + 4, gs[i].INFO.indexOf(";STD_quant_start")) - +gs[i].POS);
							if (myamt > amt) {
								mycontig = j;
								amt  = myamt;
							}
						}
					} else {
						// This segment is before the SV
						continue;
					}
				} else if ((+dgs[j].QuerySize + +dgs[j].QueryStart) < +gs[i].INFO.substring(gs[i].INFO.indexOf("END=") + 4, gs[i].INFO.indexOf(";STD_quant_start"))) {
					// This segment contains the end of the SV but not the beginning
					var myamt = (+gs[i].INFO.substring(gs[i].INFO.indexOf("END=") + 4, gs[i].INFO.indexOf(";STD_quant_start")) - +dgs[j].QueryStart) / (+gs[i].INFO.substring(gs[i].INFO.indexOf("END=") + 4, gs[i].INFO.indexOf(";STD_quant_start")) - +gs[i].POS);
					if (myamt > amt) {
						mycontig = j;
						amt  = myamt;
					}
				} else {
					// This segment is after the SV
					continue;
				}
			}
		} // end inner for
		
		// Store mycontig for drawing
		gs[i].contig = mycontig;		
	} // end outer for

	// Draw the SV rectangles
	g1.selectAll(".sv").data(wb).enter().append("rect")
		.attr("x", function(d){if (d.contig == -1) {return 0;} else {var contig = g1.selectAll("#" + dwb[+d.contig].QueryName).filter(function() {
			if (d3.select(this).attr("qstart") == dwb[d.contig].QueryStart)
				return this;
			else
				return;
		}); if (contig.size() == 0) {
				return 0;
			} else {
				return +contig.attr("x") - 2;
			}}})
		.attr("y", function(d){ if (d.contig == -1) {return 0;} else{
			var contig = g1.selectAll("#" + dwb[d.contig].QueryName).filter(function() {
				if (d3.select(this).attr("qstart") == dwb[d.contig].QueryStart)
					return this;
				else
					return;
			}); if (contig.size() == 0) {return 0;}
			var chr = contig.attr("mapchr"); 
			var endc = d.INFO.substring(d.INFO.indexOf("END=") + 4, d.INFO.indexOf(";STD_quant_start"));
			if (endc.indexOf(";") != -1) {return 0;}
			var refcoord = (((+endc - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			if (chr == "GLCHR01") {return chr1(refcoord);}
			else if (chr == "GLCHR02") {return chr2(refcoord) + c1size + chrgap;}
			else if (chr == "GLCHR03") {return chr3(refcoord);}
			else if (chr == "GLCHR04") {return chr4(refcoord) + c3size + chrgap;}
			else if (chr == "GLCHR05") {return chr5(refcoord);}
		}})
		.attr("width", 14)
		.attr("height", function(d){if (d.contig == -1) {return 0;} else {
			var contig = g1.selectAll("#" + dwb[d.contig].QueryName).filter(function() {
				if (d3.select(this).attr("qstart") == dwb[d.contig].QueryStart)
					return this;
				else
					return;
			}); if (contig.size() == 0) {return 0;}
			var chr = contig.attr("mapchr");
			var endc = d.INFO.substring(d.INFO.indexOf("END=") + 4, d.INFO.indexOf(";STD_quant_start"));
			if (endc.indexOf(";") != -1) {return 0;}
			var refcoord = (((+endc - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			var startcoord = (((+d.POS - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			if (chr == "GLCHR01") {return Math.max(1, chr1(startcoord) - chr1(refcoord));}
			else if (chr == "GLCHR02") {return Math.max(1, chr2(startcoord) - chr2(refcoord));}
			else if (chr == "GLCHR03") {return Math.max(1, chr3(startcoord) - chr3(refcoord));}
			else if (chr == "GLCHR04") {return Math.max(1, chr4(startcoord) - chr4(refcoord));}
			else if (chr == "GLCHR05") {return Math.max(1, chr5(startcoord) - chr5(refcoord));}
		}})
		.attr("class", "wb_sv")
		.attr("start", function(d){return d.POS;})
		.attr("id", function(d){return d.ID;})
		.attr("contig", function(d){return d.CHROM;})
		.attr("type", function(d){return d.ALT;})
		.attr("len", function(d){if (d.INFO.indexOf("SVLEN=") == -1){return 0;} l=d.INFO.substring(d.INFO.indexOf("SVLEN=") + 6, d.INFO.indexOf(";STRANDS")); return l;})
		.style("opacity", 0.4)
		.style("fill", "purple");
	
	g1.selectAll(".sv").data(beav).enter().append("rect")
		.attr("x", function(d){if (d.contig == -1) {return 0;} else {var contig = g1.selectAll("#" + dbeav[+d.contig].QueryName).filter(function() {
			if (d3.select(this).attr("qstart") == dbeav[d.contig].QueryStart)
				return this;
			else
				return;
		}); if (contig.size() == 0) {
				return 0;
			} else {
				return +contig.attr("x") - 2;
			}}})
		.attr("y", function(d){ if (d.contig == -1) {return 0;} else{
			var contig = g1.selectAll("#" + dbeav[d.contig].QueryName).filter(function() {
				if (d3.select(this).attr("qstart") == dbeav[d.contig].QueryStart)
					return this;
				else
					return;
			}); if (contig.size() == 0) {return 0;}
			var chr = contig.attr("mapchr"); 
			var endc = d.INFO.substring(d.INFO.indexOf("END=") + 4, d.INFO.indexOf(";STD_quant_start"));
			if (endc.indexOf(";") != -1) {return 0;}
			var refcoord = (((+endc - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			if (chr == "GLCHR01") {return chr1(refcoord);}
			else if (chr == "GLCHR02") {return chr2(refcoord) + c1size + chrgap;}
			else if (chr == "GLCHR03") {return chr3(refcoord);}
			else if (chr == "GLCHR04") {return chr4(refcoord) + c3size + chrgap;}
			else if (chr == "GLCHR05") {return chr5(refcoord);}
		}})
		.attr("width", 14)
		.attr("height", function(d){if (d.contig == -1) {return 0;} else {
			var contig = g1.selectAll("#" + dbeav[d.contig].QueryName).filter(function() {
				if (d3.select(this).attr("qstart") == dbeav[d.contig].QueryStart)
					return this;
				else
					return;
			}); if (contig.size() == 0) {return 0;}
			var chr = contig.attr("mapchr");
			var endc = d.INFO.substring(d.INFO.indexOf("END=") + 4, d.INFO.indexOf(";STD_quant_start"));
			if (endc.indexOf(";") != -1) {return 0;}
			var refcoord = (((+endc - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			var startcoord = (((+d.POS - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			if (chr == "GLCHR01") {return Math.max(1, chr1(startcoord) - chr1(refcoord));}
			else if (chr == "GLCHR02") {return Math.max(1, chr2(startcoord) - chr2(refcoord));}
			else if (chr == "GLCHR03") {return Math.max(1, chr3(startcoord) - chr3(refcoord));}
			else if (chr == "GLCHR04") {return Math.max(1, chr4(startcoord) - chr4(refcoord));}
			else if (chr == "GLCHR05") {return Math.max(1, chr5(startcoord) - chr5(refcoord));}
		}})
		.attr("class", "beav_sv")
		.attr("start", function(d){return d.POS;})
		.attr("id", function(d){return d.ID;})
		.attr("contig", function(d){return d.CHROM;})
		.attr("type", function(d){return d.ALT;})
		.attr("len", function(d){if (d.INFO.indexOf("SVLEN=") == -1){return 0;} l=d.INFO.substring(d.INFO.indexOf("SVLEN=") + 6, d.INFO.indexOf(";STRANDS")); return l;})
		.style("opacity", 0.4)
		.style("fill", "purple");
	
	g1.selectAll(".sv").data(gs).enter().append("rect")
		.attr("x", function(d){if (d.contig == -1) {return 0;} else {var contig = g1.selectAll("#" + dgs[+d.contig].QueryName).filter(function() {
			if (d3.select(this).attr("qstart") == dgs[d.contig].QueryStart)
				return this;
			else
				return;
		}); if (contig.size() == 0) {
				return 0;
			} else {
				return +contig.attr("x") - 2;
			}}})
		.attr("y", function(d){ if (d.contig == -1) {return 0;} else{
			var contig = g1.selectAll("#" + dgs[d.contig].QueryName).filter(function() {
				if (d3.select(this).attr("qstart") == dgs[d.contig].QueryStart)
					return this;
				else
					return;
			}); if (contig.size() == 0) {return 0;}
			var chr = contig.attr("mapchr"); 
			var endc = d.INFO.substring(d.INFO.indexOf("END=") + 4, d.INFO.indexOf(";STD_quant_start"));
			if (endc.indexOf(";") != -1) {return 0;}
			var refcoord = (((+endc - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			if (chr == "GLCHR01") {return chr1(refcoord);}
			else if (chr == "GLCHR02") {return chr2(refcoord) + c1size + chrgap;}
			else if (chr == "GLCHR03") {return chr3(refcoord);}
			else if (chr == "GLCHR04") {return chr4(refcoord) + c3size + chrgap;}
			else if (chr == "GLCHR05") {return chr5(refcoord);}
		}})
		.attr("width", 14)
		.attr("height", function(d){if (d.contig == -1) {return 0;} else {
			var contig = g1.selectAll("#" + dgs[d.contig].QueryName).filter(function() {
				if (d3.select(this).attr("qstart") == dgs[d.contig].QueryStart)
					return this;
				else
					return;
			}); if (contig.size() == 0) {return 0;}
			var chr = contig.attr("mapchr");
			var endc = d.INFO.substring(d.INFO.indexOf("END=") + 4, d.INFO.indexOf(";STD_quant_start"));
			if (endc.indexOf(";") != -1) {return 0;}
			var refcoord = (((+endc - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			var startcoord = (((+d.POS - +contig.attr("qstart")) * (+contig.attr("end") - +contig.attr("start"))) / (+contig.attr("qend") - +contig.attr("qstart"))) + +contig.attr("start");
			if (chr == "GLCHR01") {return Math.max(1, chr1(startcoord) - chr1(refcoord));}
			else if (chr == "GLCHR02") {return Math.max(1, chr2(startcoord) - chr2(refcoord));}
			else if (chr == "GLCHR03") {return Math.max(1, chr3(startcoord) - chr3(refcoord));}
			else if (chr == "GLCHR04") {return Math.max(1, chr4(startcoord) - chr4(refcoord));}
			else if (chr == "GLCHR05") {return Math.max(1, chr5(startcoord) - chr5(refcoord));}
		}})
		.attr("class", "gs_sv")
		.attr("start", function(d){return d.POS;})
		.attr("id", function(d){return d.ID;})
		.attr("contig", function(d){return d.CHROM;})
		.attr("type", function(d){return d.ALT;})
		.attr("len", function(d){if (d.INFO.indexOf("SVLEN=") == -1){return 0;} l=d.INFO.substring(d.INFO.indexOf("SVLEN=") + 6, d.INFO.indexOf(";STRANDS")); return l;})
		.style("opacity", 0.4)
		.style("fill", "purple");
}

/**
 * Colour things and other style elements
 */
function styleVis() {
	svg.selectAll(".wb_contig").attr("fill", "#0052f7");
	svg.selectAll(".gs_contig").attr("fill", "red");
	svg.selectAll(".beaver_contig").attr("fill", "green");
	
	// Add tooltips
	svg.selectAll(".wb_contig").on("mouseover", function(d) {
	
		var xPosition = parseFloat(d3.select(this).attr("x"));
		var yPosition = parseFloat(d3.select(this).attr("y"));

		d3.select("#contig_tooltip")
			.style("left", (xPosition + 100) + "px")
			.style("top", (yPosition + margin.top) + "px")
			.select("#contig")
			.text(d3.select(this).attr("id"));
		
		d3.select("#contig_tooltip").select("#assembly").text("AWB hybrid");
	
		d3.select("#contig_tooltip").classed("hidden", false);
	}).on("mouseout", function() {
		d3.select("#contig_tooltip").classed("hidden", true);
	});
	svg.selectAll(".beaver_contig").on("mouseover", function(d) {

		var xPosition = parseFloat(d3.select(this).attr("x"));
		var yPosition = parseFloat(d3.select(this).attr("y"));

		d3.select("#contig_tooltip")
			.style("left", (xPosition + 100) + "px")
			.style("top", (yPosition + margin.top) + "px")
			.select("#contig")
			.text(d3.select(this).attr("id"));
		
		d3.select("#contig_tooltip").select("#assembly").text("beaver hybrid");

		d3.select("#contig_tooltip").classed("hidden", false);
	}).on("mouseout", function() {
		d3.select("#contig_tooltip").classed("hidden", true);
	});
	svg.selectAll(".gs_contig").on("mouseover", function(d) {

		var xPosition = parseFloat(d3.select(this).attr("x"));
		var yPosition = parseFloat(d3.select(this).attr("y"));

		d3.select("#contig_tooltip")
			.style("left", (xPosition + 100) + "px")
			.style("top", (yPosition + margin.top) + "px")
			.select("#contig")
			.text(d3.select(this).attr("id"));
		
		d3.select("#contig_tooltip").select("#assembly").text("BGS hybrid");

		d3.select("#contig_tooltip").classed("hidden", false);
	}).on("mouseout", function() {
		d3.select("#contig_tooltip").classed("hidden", true);
	});
	
	svg.selectAll("#legend").on("mouseover", "");
	
	svg.selectAll(".wb_sv").on("mouseover", function(d) {

		var xPosition = parseFloat(+d3.select(this).attr("x") + 2);
		var yPosition = parseFloat(+d3.select(this).attr("y"));

		d3.select("#sv_tooltip")
			.style("left", (xPosition + 100) + "px")
			.style("top", (yPosition + margin.top) + "px")
			.select("#contig")
			.text(d3.select(this).attr("contig"));

		d3.select("#sv_tooltip").select("#assembly").text("AWB hybrid");
		d3.select("#sv_tooltip").select("#svid").text(d3.select(this).attr("id"));
		d3.select("#sv_tooltip").select("#type").text(d3.select(this).attr("type"));
		d3.select("#sv_tooltip").select("#len").text(d3.select(this).attr("len"));
		
		d3.select("#sv_tooltip").classed("hidden", false);
	}).on("mouseout", function() {
		d3.select("#sv_tooltip").classed("hidden", true);
	});
	svg.selectAll(".beav_sv").on("mouseover", function(d) {

		var xPosition = parseFloat(+d3.select(this).attr("x") + 2);
		var yPosition = parseFloat(+d3.select(this).attr("y"));

		d3.select("#sv_tooltip")
			.style("left", (xPosition + 100) + "px")
			.style("top", (yPosition + margin.top) + "px")
			.select("#contig")
			.text(d3.select(this).attr("contig"));

		d3.select("#sv_tooltip").select("#assembly").text("beaver hybrid");
		d3.select("#sv_tooltip").select("#svid").text(d3.select(this).attr("id"));
		d3.select("#sv_tooltip").select("#type").text(d3.select(this).attr("type"));
		d3.select("#sv_tooltip").select("#len").text(d3.select(this).attr("len"));
		
		d3.select("#sv_tooltip").classed("hidden", false);
	}).on("mouseout", function() {
		d3.select("#sv_tooltip").classed("hidden", true);
	});
	svg.selectAll(".gs_sv").on("mouseover", function(d) {

		var xPosition = parseFloat(+d3.select(this).attr("x") + 2);
		var yPosition = parseFloat(+d3.select(this).attr("y"));

		d3.select("#sv_tooltip")
			.style("left", (xPosition + 100) + "px")
			.style("top", (yPosition + margin.top) + "px")
			.select("#contig")
			.text(d3.select(this).attr("contig"));

		d3.select("#sv_tooltip").select("#assembly").text("BGS hybrid");
		d3.select("#sv_tooltip").select("#svid").text(d3.select(this).attr("id"));
		d3.select("#sv_tooltip").select("#type").text(d3.select(this).attr("type"));
		d3.select("#sv_tooltip").select("#len").text(d3.select(this).attr("len"));
		
		d3.select("#sv_tooltip").classed("hidden", false);
	}).on("mouseout", function() {
		d3.select("#sv_tooltip").classed("hidden", true);
	});

	d3.select("#selector").append("ul");
	d3.select("#selector").select("ul").selectAll("input")
		.data(sv_types).enter()
		.append("li")
		.append("label")
		.attr("class", "labels")
		.text(function(d, i){return d;})
		.append("input")
		.property("checked", true)
		.attr("type", "checkbox")
		.attr("id", function(d, i){return d;})
		.on("click", filterBoxes);
	
	d3.select("#selector").selectAll(".sizeinput").on("focusout", filterSize);

}

/**
 * Event handling for SV size filtering by user
 */
function filterSize(d, i) {
	if (d3.select(this).attr("name") == "minsize") {
		min = +d3.select(this).property("value");
		if (d3.select(this).property("value") == "") {
			min = 0;
		}
	} else if (d3.select(this).attr("name") == "maxsize") {
		max = +d3.select(this).property("value");
		if (d3.select(this).property("value") == "") {
			max = 4471437;
		}
	}
	
	svg.selectAll(".wb_sv").filter(function() {
		if(d3.select(this).attr("type_hidden")) {
			return 0;
		} else {
			return 1;
		}
	}).attr("width", 14).attr("size_hidden", false);
	svg.selectAll(".beav_sv").filter(function() {
		if(d3.select(this).attr("type_hidden")) {
			return 0;
		} else {
			return 1;
		}
	}).attr("width", 14).attr("size_hidden", false);
	svg.selectAll(".gs_sv").filter(function() {
		if(d3.select(this).attr("type_hidden")) {
			return 0;
		} else {
			return 1;
		}
	}).attr("width", 14).attr("size_hidden", false);
	
	svg.selectAll(".wb_sv").filter(function() {
		if(Math.abs(+d3.select(this).attr("len")) < min || Math.abs(+d3.select(this).attr("len")) > max) {
			return 1;
		} else {
			return 0;
		}
	}).attr("width", "0").attr("size_hidden", true);
	svg.selectAll(".beav_sv").filter(function() {
		if(Math.abs(+d3.select(this).attr("len")) < min || Math.abs(+d3.select(this).attr("len")) > max) {
			return 1;
		} else {
			return 0;
		}
	}).attr("width", "0").attr("size_hidden", true);
	svg.selectAll(".gs_sv").filter(function() {
		if(Math.abs(+d3.select(this).attr("len")) < min || Math.abs(+d3.select(this).attr("len")) > max) {
			return 1;
		} else {
			return 0;
		}
	}).attr("width", "0").attr("size_hidden", true);
}

/**
 * Event handling function for SV filtering by user
 */
function filterBoxes(d, i) {
	var type = "";
	
	switch(d) {
	case "Duplications":
		type = "<DUP>";
		break;
	case "Deletions":
		type = "<DEL>";
		break;
	case "Insertions":
		type = "<INS>";
		break;
	case "Inversions":
		type = "<INV>";
		break;
	case "Inverted Duplications":
		type = "<INVDUP>";
		break;
	}

	d3.selectAll("input").filter(function(){
		if(d3.select(this).attr("id") == d) {
			if(d3.select(this).property("checked") == true)
				svg.selectAll("rect").filter(function() {
					if(d3.select(this).attr("size_hidden")) {
						return 0;
					} else {
						return 1;
					}
				}).filter(function() {
					if (d3.select(this).attr("type") == type)
						return 1;
					else 
						return 0;
				}).attr("width", 14).attr("type_hidden", false);
			else
				svg.selectAll("rect").filter(function() {
					if(d3.select(this).attr("size_hidden")) {
						return 0;
					} else {
						return 1;
					}
				}).filter(function() {
					if (d3.select(this).attr("type") == type)
						return 1;
					else 
						return 0;
				}).attr("width", 0).attr("type_hidden", true);
			
		}
	});
}
