// This recieves messages of type "setActiveTab"
Shiny.addCustomMessageHandler("setActiveTab",
	function(message) {
		console.log("called");
		console.log(message);
		tabs = $('.nav li')
		tabs.each(function() {
			$(this).removeClass('active')
		})

		// forces cluster, this needs to be fixed
		$(tabs[message.tabNo]).addClass('active')
		tabsContents = $('.tab-content .tab-pane')
		tabsContents.each(function() {
			$(this).removeClass('active')
		})
		// forces cluster, this needs to be fixed
		$(tabsContents[message.tabNo]).addClass('active')
		// forces show of cluster drop down, this needs to be fixed
		$(message.tabControl).trigger('change').trigger('shown');
	}
);

var restore_snapshot = function(el, input_params) {
		/* Shiny.inputBindings.getBindings() return the InputBinding instances
		   for every (native) input type that Shiny supports (selectInput, textInput,
		   actionButton etc.)  */
		$.each(Shiny.inputBindings.getBindings(), function(i, b) {
			/* find all inputs within a specific input type */
			var inputs = b.binding.find(el);
			$.each(inputs, function(j, inp) {
				/* check if the input's id matches the key specified in the query
				   string */
				var inp_val = input_params[$(inp).attr("id")];
				console.log($(inp).attr("id"));
				if (inp_val != undefined) {
					console.log(inp_val);
					b.binding.setValue(inp, inp_val);
				}
			});
		});
		// force some updates
		$("#myClusterGenesUpdateButton").click();
	}

var restoreBinding = new Shiny.OutputBinding();
	$.extend(restoreBinding, {
		find: function(scope) {
			return $(scope).find("#inputContainer");
		},
		renderValue: function(el, data) {
			console.log("data:");
			console.log(data);
			console.log(":data");
			// very rudimentary sanity check
			if ($.isPlainObject(data) && data.hasOwnProperty('k')) {
				restore_snapshot(el, data);
				alert("Snapshot restored!");
			}
		}
	});

 Shiny.outputBindings.register(restoreBinding, 'inputs.Restore');
