function flag = hasMenu(obj, menuKey)
    %HASMENU Return true if we have a menu item by key
    flag = ischar(menuKey) && isKey(obj.hMenus, menuKey);
end