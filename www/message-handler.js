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
