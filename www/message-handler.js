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

Shiny.addCustomMessageHandler("myClusterGenesUpdateButtonClick",
	function(message) {
		$("#myClusterGenesUpdateButton").click();
	}
);

Shiny.addCustomMessageHandler("testMsg",
	function(message) {
		console.log("hai!")
		$("#k").value = 40;
		$("#cluster").value = 12;
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
                if (inp_val != undefined) {
                    b.binding.setValue(inp, inp_val);
                }
            });
        });
    }

var restoreBinding = new Shiny.OutputBinding();
    $.extend(restoreBinding, {
        find: function(scope) {
            return $(scope).find("#input_container");
        },
        renderValue: function(el, data) {
            // very rudimentary sanity check
			console.log(data)
            if ($.isPlainObject(data) && data.hasOwnProperty('username')) {
                restore_snapshot(el, data);
                alert("Snapshot restored!");
            }
        }
    });

    Shiny.outputBindings.register(restoreBinding, 'inputs.Restore');
