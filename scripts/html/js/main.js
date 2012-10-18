$(document).ready(function(){

	// Where am I?
	var mypath = location.href;
	var str = mypath.split("/");
	var estr = str[0];
	for (i=1;i<str.length-1;i++) {
		estr = estr + "/" + str[i];
	}

	// Get JSON
	var chromjson = estr+"/conf/chrom.json";	
	var samplejson = estr+"/conf/samples.json";

        $.ajaxSetup ({ cache: false});

	jQuery.support.cors = true;

	$.getJSON(chromjson, function(data) {
		$.each(data, function(key, val) {
			var listin1 = $("#chrom1");
			var listin2 = $("#chrom2");
			$('<option value="'+val+'">'+val+'</option>').prependTo(listin1);
			$('<option value="'+val+'">'+val+'</option>').prependTo(listin2);
		});
	});
	
	
	$.getJSON(samplejson, function(data) {

		$.each(data, function(key, val) {
			var samplesin = $("#samples");

			$('<option value="'+val['id']+'">'+val['id']+'</option>').prependTo(samplesin);
		});
	});
	
	$("#samples").change(function(){
		
		var samplevalue = $(this).val();
		
		$('#sets').find('option').remove().end().append('<option value="" selected="selected">Select a sample</option>').val('');

		$.getJSON(samplejson, function(data) {

			$.each(data, function(key, val) {
				var setsin = $("#sets");
				
				if (val['id'] == samplevalue) {
					
					$.each(val['value'], function(i, o) {
						$('<option value="'+o+'">'+o+'</option>').prependTo(setsin);

					});
				}

			});
		});
	});
 
        $("#chrom1").change(function(){
		
                        var value = $(this).val();
                        
                        var value2 = $("#chrom2").val();
			
			if (value == 'all' && value2 == 'all') {
				
				//Do different stuff; show all options
				makevisible($("#listexps"));
				makevisible($(".auto0"));
				makevisible($(".auto1"));
			}
                        
			else {
			
				if (value!='' && value2!='') {
				    makevisible($("#listexps"));
				    
				    if (value != value2) {
					makevisible($(".auto0"));
					makehidden($(".auto1"));
				    }
				    
				    else {
					makevisible($(".auto1"));
					makehidden($(".auto0"));
				    }
				}
				
			}

        });

        $("#chrom2").change(function(){
		
                        var value2 = $(this).val();
                        
                        var value = $("#chrom1").val();
                       
			if (value == 'all' && value2 == 'all') {
				
				//Do different stuff; show all options
				makevisible($("#listexps"));
				makevisible($(".auto0"));
				makevisible($(".auto1"));
			}
			
			else {
			
				if (value!='' && value2!='') {
				    
				    makevisible($("#listexps"));
				    
				    if (value != value2) {
					makevisible($(".auto0"));
					makehidden($(".auto1"));
				    }
				    
				    else {
					makevisible($(".auto1"));
					makehidden($(".auto0"));
				    }
				}
			}
        });
        
	


         
        $("a.exp").click(function(){
    
           var experiment = $(this).text();
           var chrom1 = $("#chrom1").val();
           var chrom2 = $("#chrom2").val();
	   var setexp = $("#sets").val();
           var setsam = $("#samples").val();

	   
            if (chrom1.length>0 && chrom2.length>0 && setexp.length>0) {
        
                var resultat = $("#result");
                addboxes(resultat, experiment, chrom1, chrom2, setexp, setsam);
		$('.warning').remove();
		
            
	    }
            
	    else {
			var warning = $("#container");
			$('<p class="warning">SELECT ALL THE OPTIONS</p>').prependTo(warning);
		
	    }
	    
            return(false);
        });




});


function clearelement(element) {
    
    $("#"+element).remove();
}

function clearelementall(element1, element2) {
    
    $("."+element1).remove();
    $("."+element2).remove();
    $('.clearall').hide();
}

function makevisible (element) {

    element.removeClass("hidden");

}

function makehidden (element) {
    
    element.addClass("hidden");

}



function addboxes (result, experiment, chrom1, chrom2, setexp, setsam) {
    var randomNum = Math.ceil(Math.random()*1000);
    
    //Display clear all
    $('.clearall').show();
    
    // Where am I?
    var mypath = location.href;
    var str = mypath.split("/");
    var estr =str[0];
    for (i=1;i<str.length-1;i++) {
	estr = estr +"/"+ str[i];
    }
	
    // Get JSON
    var chromjson = estr+"/conf/chrom.json";
    

    $.getJSON(chromjson, function(data) {

     var items = [];

	$.each(data, function(key, val) {
	items.push(val);
	});

    var middle = setexp.replace('bin.', '');
    var dirresults = "results."+setsam;
    
    var intercorresp = new Array();
    intercorresp['correlation'] = 'heatmap.correlation.png';
    intercorresp['coverage-correction'] = 'heatmap.cvg.png';
    intercorresp['observed'] = 'heatmap.observed.png';
    intercorresp['expected'] = 'heatmap.expected.png';
    intercorresp['observed.vs.expected'] = 'heatmap.intra-interaction.png';

    if ((chrom1 == 'all') || (chrom2 == 'all')) {
	
		if ( chrom1 == chrom2 ) {
			
			$('.warning').remove(); //Remove any msg
			var warning = $("#container");
			$('<p class="warning">Selecting ALL CHROMOSOMES will take a while!</p>').prependTo(warning);
			
			$.each(items, function(key1, value1) {
				var chromvalue1 = value1;
				$.each(items, function(key2, value2) {
					var chromvalue2 = value2;
					//cases of no repeat
					if (experiment == 'expected' || experiment == 'observed.vs.expected') {
						if (chromvalue2 != chromvalue1) {
							//Skip if not equal
							return true;
						}
					}
					$('<div id="output-'+randomNum+'" class="output-'+experiment+' outputlittle"><p class="clearit"><a name="a-'+experiment+'" /><a href="#a-'+experiment+'" onclick="clearelement(\'output-'+randomNum+'\')">Clear</a></p><p class="output-title">'+chromvalue1+'-'+chromvalue2+' ('+experiment+')</p><p class="output-title2">'+setsam+' - '+setexp+'</p><a href="'+dirresults+'/'+setexp+'/'+experiment+'/'+chromvalue1+'-'+chromvalue2+'.'+middle+'.'+intercorresp[experiment]+'"><img src="'+dirresults+'/'+setexp+'/'+experiment+'/'+chromvalue1+'-'+chromvalue2+'.'+middle+'.'+intercorresp[experiment]+'" /></a></div>').prependTo(result);
				});
			});
			
		}
		
		else {
			
			$('.warning').remove(); //Remove any msg
			//Get all values from chromosome
			if (chrom1 == 'all') {
				$.each(items, function(key, value) {
					var chromvalue = value;
					$('<div id="output-'+randomNum+'" class="output-'+experiment+' output"><p class="clearit"><a name="a-'+experiment+'" /><a href="#a-'+experiment+'" onclick="clearelement(\'output-'+randomNum+'\')">Clear</a></p><p class="output-title">'+chromvalue+'-'+chrom2+' ('+experiment+')</p><p class="output-title2">'+setsam+' - '+setexp+'</p><a href="'+dirresults+'/'+setexp+'/'+experiment+'/'+chromvalue+'-'+chrom2+'.'+middle+'.'+intercorresp[experiment]+'"><img src="'+dirresults+'/'+setexp+'/'+experiment+'/'+chromvalue+'-'+chrom2+'.'+middle+'.'+intercorresp[experiment]+'" /></a></div>').prependTo(result);

				});
			}
			
			//Get all values from chromosome
			if (chrom2 == 'all') {
				$.each(items, function(key, value) {
					var chromvalue = value;
					$('<div id="output-'+randomNum+'" class="output-'+experiment+' output"><p class="clearit"><a name="a-'+experiment+'" /><a href="#a-'+experiment+'" onclick="clearelement(\'output-'+randomNum+'\')">Clear</a></p><p class="output-title">'+chrom1+'-'+chromvalue+' ('+experiment+')</p><p class="output-title2">'+setsam+' - '+setexp+'</p><a href="'+dirresults+'/'+setexp+'/'+experiment+'/'+chrom1+'-'+chromvalue+'.'+middle+'.'+intercorresp[experiment]+'"><img src="'+dirresults+'/'+setexp+'/'+experiment+'/'+chrom1+'-'+chromvalue+'.'+middle+'.'+intercorresp[experiment]+'" /></a></div>').prependTo(result);

				});
			}
			
		}
    }
    
    else {
    
	$('.warning').remove(); //Remove any msg
	$('<div id="output-'+randomNum+'" class="output-'+experiment+' output"><p class="clearit"><a name="a-'+experiment+'" /><a href="#a-'+experiment+'" onclick="clearelement(\'output-'+randomNum+'\')">Clear</a></p><p class="output-title">'+chrom1+'-'+chrom2+' ('+experiment+')</p><p class="output-title2">'+setsam+' - '+setexp+'</p><a rel="chromimg" href="'+dirresults+'/'+setexp+'/'+experiment+'/'+chrom1+'-'+chrom2+'.'+middle+'.'+intercorresp[experiment]+'"><img src="'+dirresults+'/'+setexp+'/'+experiment+'/'+chrom1+'-'+chrom2+'.'+middle+'.'+intercorresp[experiment]+'" /></a></div>').prependTo(result);
    }

	});
	return(true);
}


