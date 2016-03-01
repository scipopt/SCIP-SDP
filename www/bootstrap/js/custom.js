! function ($) {
    $(function () {
        var $window = $(window);
        var $body = $(document.body);
        var $sideBar = $('.bs-sidebar');
        var navHeight = $('.navbar').outerHeight(true) + 10;

        $(".answer").hide();
        var elem = window.location.hash.replace('#', '');
        if(elem) {
            // scroll to correct position
            location.href = "#"+elem;
            // reveal answer if required
            target = $("#"+elem);
            if (target.hasClass("reveal_faq")) {
                var answer = $(href+"_ans");
                answer.show();
            }
        }

        // reveal faq answer
        $(".reveal_faq").click(function(event) {
            // construct answer id
            var questionId = "#" + $(event.target).parent().prop("id") + "_ans";
            // reveal style
            $(questionId).toggle("fast");
        });

        $(".reveal").click(function(event) {
            // construct answer id
            var questionId = "#" + $(event.target).prop("id") + "_ans";
            // reveal style
            $(questionId).toggle("fast");
        });

        // reveal linked answer
        $(".link-answer").click(function(event) {
            var href = $ (event.target).attr('href');
            target = $(href);
            if ((target.hasClass("reveal")) || (target.hasClass("reveal_faq"))) {
                var answer = $(href+ "_ans");
                answer.show();
            }
        });

        $body.scrollspy({
            target: '.bs-sidebar',
            offset: navHeight
        });

        $('.bs-docs-container [href=#]').click(function (e) {
            e.preventDefault()
        });

        $window.on('resize', function () {
            $body.scrollspy('refresh')
            // We were resized. Check the position of the nav box
            $sideBar.affix('checkPosition')
        });

        $window.on('load', function () {
            $body.scrollspy('refresh')
            $('.bs-top').affix();
            $sideBar.affix({
                offset: {
                    top: function () {
                        var offsetTop = $sideBar.offset().top
                        var sideBarMargin = parseInt($sideBar.children(0).css('margin-top'), 10)
                        var navOuterHeight = $('.bs-docs-nav').height()

                        // We can cache the height of the header (hence the this.top=)
                        // This function will never be called again.
                        return (this.top = offsetTop - navOuterHeight - sideBarMargin);
                    },
                    bottom: function () {
                        // We can't cache the height of the footer, since it could change
                        // when the window is resized. This function will be called every
                        // time the window is scrolled or resized
                        return $('.bs-footer').outerHeight(true)
                    }
                }
            })
            setTimeout(function () {
                // Check the position of the nav box ASAP
                $sideBar.affix('checkPosition')
            }, 10);
            setTimeout(function () {
                // Check it again after a while (required for IE)
                $sideBar.affix('checkPosition')
            }, 100);
        });

        // tooltip demo
        $('.tooltip-demo').tooltip({
            selector: "[data-toggle=tooltip]",
            container: "body"
        });

        $('.tooltip-test').tooltip();
        $('.popover-test').popover();

        $('.bs-docs-navbar').tooltip({
            selector: "a[data-toggle=tooltip]",
            container: ".bs-docs-navbar .nav"
        });
    })
}(window.jQuery)