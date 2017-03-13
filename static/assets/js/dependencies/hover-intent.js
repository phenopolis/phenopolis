/*!
* hoverIntent v1.8.1 // 2014.08.11 // jQuery v1.9.1+
* http://briancherne.github.io/jquery-hoverIntent/
*
* You may use hoverIntent under the terms of the MIT license. Basically that
* means you are free to use hoverIntent as long as this header is left intact.
* Copyright 2007, 2014 Brian Cherne
*/

!function(e){"use strict";"function"==typeof define&&define.amd?define(["jquery"],e):jQuery&&!jQuery.fn.hoverIntent&&e(jQuery)}(function(e){"use strict";var t,n,o={interval:100,sensitivity:6,timeout:0},i=0,u=function(e){t=e.pageX,n=e.pageY},r=function(e,o,i,v){return Math.sqrt((i.pX-t)*(i.pX-t)+(i.pY-n)*(i.pY-n))<v.sensitivity?(o.off(i.event,u),delete i.timeoutId,i.isActive=!0,e.pageX=t,e.pageY=n,delete i.pX,delete i.pY,v.over.apply(o[0],[e])):(i.pX=t,i.pY=n,i.timeoutId=setTimeout(function(){r(e,o,i,v)},v.interval),void 0)},v=function(e,t,n,o){return delete t.data("hoverIntent")[n.id],o.apply(t[0],[e])};e.fn.hoverIntent=function(t,n,a){var s=i++,d=e.extend({},o);e.isPlainObject(t)?(d=e.extend(d,t),e.isFunction(d.out)||(d.out=d.over)):d=e.isFunction(n)?e.extend(d,{over:t,out:n,selector:a}):e.extend(d,{over:t,out:t,selector:n});var f=function(t){var n=e.extend({},t),o=e(this),i=o.data("hoverIntent");i||o.data("hoverIntent",i={});var a=i[s];a||(i[s]=a={id:s}),a.timeoutId&&(a.timeoutId=clearTimeout(a.timeoutId));var f=a.event="mousemove.hoverIntent.hoverIntent"+s;if("mouseenter"===t.type){if(a.isActive)return;a.pX=n.pageX,a.pY=n.pageY,o.off(f,u).on(f,u),a.timeoutId=setTimeout(function(){r(n,o,a,d)},d.interval)}else{if(!a.isActive)return;o.off(f,u),a.timeoutId=setTimeout(function(){v(n,o,a,d.out)},d.timeout)}};return this.on({"mouseenter.hoverIntent":f,"mouseleave.hoverIntent":f},d.selector)}});
