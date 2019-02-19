function preview(obj)
    %PREVIEW Display Preview GUI
    obj.startParPool();
    hPreview = jrclust.views.PreviewController(obj.hCfg);
    hPreview.preview();
end