// This recieves messages of type "restore.like.state" from the server.
Shiny.addCustomMessageHandler("restore.like.state",
	function(message) {
		for (property in message) {
			var obj = document.getElementById(property);
			if (obj != null) {
				Shiny.onInputChange(property, message[property][0]);
				//obj.value = message[property][0];
				//console.log(message[property][0]);
				//console.log(obj.value);
				//console.log(property);
				//console.log(obj);
			}
		}
		//console.log(message);
		//console.log(JSON.stringify(message));
		// this should now switch the tab
		tabs = $('.nav li')
		tabs.each(function() {
			$(this).removeClass('active')
		})

		// forces cluster, this needs to be fixed
		$(tabs[2]).addClass('active')
		tabsContents = $('.tab-content .tab-pane')
		tabsContents.each(function() {
			$(this).removeClass('active')
		})
		// forces cluster, this needs to be fixed
		$(tabsContents[2]).addClass('active')
		// forces show of cluster drop down, this needs to be fixed
		$('#cluster').trigger('change').trigger('shown');
	}
);
