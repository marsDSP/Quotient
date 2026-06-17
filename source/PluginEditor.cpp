#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PluginProcessorEditor::PluginProcessorEditor (PluginProcessor& p)
    : AudioProcessorEditor (&p), pRef (p)
{
    ignoreUnused (pRef);
    setSize (400, 300);
}

PluginProcessorEditor::~PluginProcessorEditor() = default;

//==============================================================================
void PluginProcessorEditor::paint (Graphics& g)
{
}

void PluginProcessorEditor::resized()
{
}
